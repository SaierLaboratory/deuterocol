#!/usr/bin/env python2

from __future__ import print_function, division, generators

import argparse, json, os, subprocess, re, sys, time
import tempfile
import superpose
import shutil

import tmalignparser


VERBOSITY = 1
def info(*things):
	print('[INFO]:', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]:', *things, file=sys.stderr)
def error(*things):
	print('[ERROR]:', *things, file=sys.stderr)
	exit(1)

def count_residues(pdbfn):
	n = 0
	with open(pdbfn) as f:
		for l in f:
			if l.startswith('ATOM'): 
				if l[13:15] == 'CA': n += 1
	return n

def _intersection_size(spans1, spans2):
	n = 0
	lastn = None
	for span1 in spans1:
		for span2 in spans2:
			if span1[0] <= span2[0] <= span2[1] <= span1[1]: 
				n += span2[1] - span2[0] + 1
			elif span2[0] <= span1[0] <= span1[1] <= span2[1]: 
				n += span1[1] - span1[0] + 1
			elif (span2[0] == span1[0]) and (span1[1] == span2[1]): 
				n += span1[1] - span1[0] + 1
			elif span1[0] <= span2[0] <= span1[1] or span2[0] <= span1[1] <= span2[1]: 
				n += span1[1] - span2[0] + 1
			elif span2[0] <= span1[0] <= span2[1] or span1[0] <= span2[1] <= span1[1]:
				n += span2[1] - span1[0] + 1
			else: 
				n += 0
			if lastn is not None and n < lastn: 
				print(span1, span2)
		lastn = n
	return n

class TMalign(superpose.Superpose):
	def __init__(self, d2dir='deuterocol2', force=False, skip_cut=True):
		superpose.Superpose.__init__(self, d2dir=d2dir, force=force)
		self.skip_cut = skip_cut

	def main(self, famdir):
		done = []
		if not os.path.isdir('{}/tmalignments'.format(famdir)): os.mkdir('{}/tmalignments'.format(famdir))
		if VERBOSITY: info('Checking for existing alignments in {}...'.format(famdir))
		if not self.force and os.path.isfile('{}/tmalignments/sp_all.tsv'.format(famdir)):
			with open('{}/tmalignments/sp_all.tsv'.format(famdir)) as f:
				for l in f:
					if not l.strip(): continue
					elif l.startswith('#'): continue
					else: 
						try: 
							json.loads(l.split('\t')[1])
							done.append(l.split('\t')[0])
						except ValueError: break
		if done: info('Skipping {} alignments (already done)'.format(len(done)))

		todo = -len(done)
		with open('{}/config/agenda.json'.format(famdir)) as g:
			for l in g: 
				if not l.strip(): continue
				elif l.startswith('#'): continue
				todo += 1
		n = 0
		t = superpose.time.time()
		if VERBOSITY: info('Performing {} alignments with TM-align...'.format(todo))
		if todo:
			with open('{}/tmalignments/sp_all.tsv'.format(famdir), 'a') as f:
				with open('{}/config/agenda.json'.format(famdir)) as g:
					for l in g:
						alnstart = time.time()
						obj = json.loads(l)
						if not self.force and obj['name'] in done: 
							n += 1
							continue
						if not n % 100:
							info('Finished {}/{} alignments({:0.2f}s since last message)'.format(n, todo, superpose.time.time() - t))
							t = superpose.time.time()
						query, subject = obj['name'].split('_vs_')
						tf = tempfile.NamedTemporaryFile()
						qfn = '{}/../cut_pdbs/{}.pdb'.format(famdir, query)
						sfn = '{}/../cut_pdbs/{}.pdb'.format(famdir, subject)
						cmd = ['TMalign', 
							qfn, sfn, 
							'-m', tf.name
						]

						#this is where TMalign is run
						p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
						out, err = p.communicate()
						if p.returncode:
							#29: file not found, e.g. from being incompetently eliminated from the pipeline
							#174: segfault, seems to be from having empty structures
							if p.returncode in {29, 174}: 
								n += 1
								continue
							else: raise Exception(err)
						tf.flush()
						tf.seek(0)

						#capture stderr to do this later #TODO
						#if out.startswith('forrtl: No such file'): 
						#	superpose.error('Could not find at least one of {} and {}'.format(qfn, tfn))
						#for l in out.split('\n'): print(l)
						sp = superpose.Alignment()
						for i, rawrow in enumerate(tf):
							if 2 <= i <= 4:
								row = [float(x) for x in rawrow.strip().split()[1:]]
								sp.matrix.append(row[1:] + row[:1])
						sp.qpresent = [obj['qspan']]
						sp.spresent = [obj['sspan']]
						sp.query = obj['query']
						sp.subject = obj['subject']
						sp.queryfn = qfn
						sp.subjectfn = sfn
						#print(sp.dump_json())
						tmscores = []
						alignment = False
						aligned = []

						#this is the where TMalign parsing starts
						for tml in out.split('\n'):
							if alignment is False:
								if tml.startswith('Aligned length'):
									sp.length = int(tml.split()[2][:-1])
									sp.rmsd = float(tml.split()[4][:-1])
								elif tml.startswith('TM-score'):
									tmscores.append(float(tml.split()[1]))
								elif tml.startswith('(":" denotes'): alignment = True
							else:
								if tml.strip():
									aligned.append(tml.replace('\n', ''))
						sp.quality = max(tmscores)
						qseq, dseq, sseq = aligned

						covstuff = True
						if covstuff:
							covstart = time.time()

							#if a contig couldn't be found but the fragment was shorter than this, skip
							minl_relevant = 4

							qatomseq = tmalignparser.extract_atom_sequence(open(qfn), obj['qchain'])
							satomseq = tmalignparser.extract_atom_sequence(open(sfn), obj['schain'])

							permitted = '.:'
							lastd = '#'
							qfrags = []
							sfrags = []
							for q, d, s in zip(qseq, dseq, sseq):
								if d in permitted: 
									if lastd in permitted:
										qfrags[-1] += q
										sfrags[-1] += s
									else:
										qfrags.append(q)
										sfrags.append(s)
								lastd = d

							qistart = None
							qaln_combined = tmalignparser.NumberedSequence()
							for qfrag in qfrags:
								qcontig = qatomseq.gapless_align_string(qfrag, start=qistart)
								if qcontig is None: 
									if len(qfrag) < minl_relevant: continue
									#print(qfrag, qatomseq)
									continue

								qaln_combined = qaln_combined + qcontig
								qistart = qcontig.get_range()[-1] + 1

							sistart = None
							saln_combined = tmalignparser.NumberedSequence()
							for sfrag in sfrags:
								scontig = satomseq.gapless_align_string(sfrag, start=sistart)
								#print(sfrag, sistart, scontig, satomseq)
								if scontig is None: 
									if len(sfrag) < minl_relevant: continue
									#print(sfn)
									#print(satomseq.iterable)
									#print(out)
									#exit(1)
									continue

								saln_combined = saln_combined + scontig
								#sistart = scontig.get_range()[-1] + 1
								sistart = scontig.get_range()[-1] - 2

							n_qaligned = _intersection_size(qaln_combined.get_ranges(), obj['qindices'])
							n_saligned = _intersection_size(saln_combined.get_ranges(), obj['sindices'])

							sp.qlen = len(qatomseq)
							sp.slen = len(satomseq)

							sp.qtmlen = obj['qlen']
							sp.stmlen = obj['slen']

							sp.qtmcov = n_qaligned / sp.qtmlen
							sp.stmcov = n_saligned / sp.stmlen
							sp.tmcov = max(sp.qtmcov, sp.stmcov)

							sp.qfullcov = sp.length / sp.qlen
							sp.sfullcov = sp.length / sp.slen
							sp.fullcov = max(sp.qfullcov, sp.sfullcov)


							covtime = time.time() - covstart
						dseq += ' ' * max(0, len(qseq)-len(dseq))
						qaligned, saligned, distances = [], [], []
						aligned = 0


						lastqr, lastsr = '', ''
						lastdist = -1
						qi, si = sp.qpresent[0][0], sp.spresent[0][0]
						for qcur, dist, scur in zip(qseq, dseq, sseq):
							if lastdist  == -1:
								if dist == ' ': 
									qaligned.append([None, 1])
									saligned.append([None, 1])
									distances.append(None)
								elif dist == '.':
									if qcur == '-': qaligned.append([None, 1])
									else: qaligned.append([qi, 1])
									if scur == '-': saligned.append([None, 1])
									else: saligned.append([si, 1])
									distances.append(5.0)
									aligned += 1
								elif dist == ':':
									if qcur == '-': qaligned.append([None, 1])
									else: qaligned.append([qi, 1])
									if scur == '-': saligned.append([None, 1])
									else: saligned.append([si, 1])
									distances.append(0.0)
									aligned += 1
							elif dist == ' ': 
								if qaligned[-1][0] is None: qaligned[-1][1] += 1
								else: qaligned.append([None, 1])
								if saligned[-1][0] is None: saligned[-1][1] += 1
								else: saligned.append([None, 1])
								distances.append(None)
							elif dist == '.': 
								if qaligned[-1][0] is None: qaligned.append([qi, 1])
								else: qaligned[-1][1] += 1
								if saligned[-1][0] is None: saligned.append([si, 1])
								else: saligned[-1][1] += 1
								distances.append(5.0)
								aligned += 1
							elif dist == ':': 
								if qaligned[-1][0] is None: qaligned.append([qi, 1])
								else: qaligned[-1][1] += 1
								if saligned[-1][0] is None: saligned.append([si, 1])
								else: saligned[-1][1] += 1
								distances.append(0.)
								aligned += 1
							#if not qaligned and qcur == '-': qaligned.append([None, 1])
							#elif not qaligned and qcur != '-': qaligned.append([qi, 1])
							#elif qcur == '-' and lastqr == '-': qaligned[-1][1] += 1
							#elif qcur == '-' and lastqr != '-': qaligned.append([None, 1])
							#elif qcur != '-' and lastqr != '-': qaligned[-1][1] += 1
							#elif qcur != '-' and lastqr == '-': qaligned.append([qi, 1])

							#if not saligned and scur == '-': saligned.append([None, 1])
							#elif not saligned and scur != '-': saligned.append([si, 1])
							#elif scur == '-' and lastsr == '-': saligned[-1][1] += 1
							#elif scur == '-' and lastsr != '-': saligned.append([None, 1])
							#elif scur != '-' and lastsr != '-': saligned[-1][1] += 1
							#elif scur != '-' and lastsr == '-': saligned.append([si, 1])
						
							##TODO: expose distance cutoff
							#if dist == ' ': distances.append(None)
							#elif dist == '.': distances.append(5.0)
							#elif dist == ':': distances.append(0.)

							if qcur != '-': qi += 1
							if scur != '-': si += 1
							if qcur == '-': qi += 1
							if scur == '-': si += 1
							lastqr = qcur
							lastsr = scur
							lastdist = dist

						sp.qaligned = qaligned
						sp.saligned = saligned
						sp.distances = distances
						sp.aligned = aligned

						query, qchain, qhel, vs, subject, schain, shel = obj['name'].split('_')
						sp.qhel = [int(x) for x in qhel[1:].split('-')]
						sp.shel = [int(x) for x in shel[1:].split('-')]
						f.write('{}\t{}\n'.format(obj['name'], sp.dump_json()))
						n += 1
						alntime = time.time() - alnstart
						#print(covtime, alntime)
		info('Finished alignments for {}'.format(famdir))


	def cut_pdbs(self, famdir):
		indices = {}
		with open('{}/config/indices.json'.format(famdir)) as f:
			for l in f:
				sl = l.split('\t')
				indices[sl[0]] = json.loads(sl[1])

		bundle = -1
		with open('{}/config/agenda.json'.format(famdir)) as f:
			obj = json.loads(f.readline())
			bundle = obj['bundle']

		cutme = set()
		with open('{}/config/align_me.json'.format(famdir)) as f:
			for l in f:
				obj = json.loads(l)
				for fam in obj:
					for pdbid in obj[fam]:
						cutme.add(pdbid[:4])

		chains = {}
		for pdbc in indices:
			try: chains[pdbc[:4]].append(pdbc)
			except KeyError: chains[pdbc[:4]] = [pdbc]
		
		for pdb in sorted(cutme):
			for pdbc in chains[pdb]:
				if len(indices[pdbc]) <= 2: continue
				#if bundle >= len(indices[pdbc]): 
				if False:
					shutil.copy('{}/../pdbs/{}.pdb'.format(famdir, pdb), 
						'{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdbc, 1, len(indices[pdbc])))
				else:
					for start in range(0, len(indices[pdbc])-bundle+1):
						end = start + bundle - 1
						infile = '{}/../pdbs/{}.pdb'.format(famdir, pdb)
						outfile = '{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdbc, start+1, end+1)
						if not os.path.isfile(infile): continue
						if os.path.isfile(outfile) and os.path.getsize(outfile):
							with open(outfile) as f:
								for l in f: pass
							if l.startswith('END'): continue


						satisfied = False
						#needed because endpoints found in unsolved residues prevent cutting in their direction
						theorleft = indices[pdbc][end][1] - indices[pdbc][start][0] + 1
						i = 0
						headcorrection = 0
						tailcorrection = 0
						while not satisfied:
							i += 1
							#this version requires manual inspection for residue counts but somehow doesn't segfault
							cmd = 'lvresidue {}-{}\nlvchain {}\nwrite PDB'.format(
								indices[pdbc][start][0]+headcorrection,
								indices[pdbc][end][1]-tailcorrection,
								pdbc[-1])

							#this version segfaults for some reason
							#cmd = 'lvchain {}\nlvresidue {}-{}\nwrite PDB'.format(
							#	pdbc[-1],
							#	indices[pdbc][start][0]+headcorrection,
							#	indices[pdbc][end][1]-tailcorrection,
							#	)
							p = subprocess.Popen(['pdbcur', 
								'xyzin', infile,
								'xyzout', outfile],
								stdin=subprocess.PIPE,
								stdout=subprocess.PIPE,
								stderr=subprocess.PIPE)
							out, err = p.communicate(input=cmd)
							#if VERBOSITY and out.strip(): print(out)

							try: empleft = count_residues(outfile)
							except IOError: empleft = 0
							#print(pdbc, empleft, theorleft, empleft > (theorleft + 5))

							if empleft > (theorleft + 5):
								if i % 2: tailcorrection += 1
								#else: headcorrection += 1
								else: tailcorrection += 1
							else: satisfied = True
						if headcorrection or tailcorrection:
							warn('{} for helices {}-{} has TM indices assigned on gaps. Chopped off {} residues to make the bug go away'.format(pdbc, start+1, end+1, headcorrection+tailcorrection))
							#print('{} {}-{} ({}/{}, chopped N-{} C-{})'.format(pdbc, indices[pdbc][start][0], indices[pdbc][end][1], empleft, theorleft, headcorrection, tailcorrection))


	def run(self):
		if VERBOSITY: info('Checking for expected directory contents...')
		for famdir in self.famdirs:
			if not os.path.isfile('{}/config/agenda.json'.format(famdir)): raise IOError('Missing agenda.json: Deuterocol2 directory structure not complete for {}'.format(famdir))
			if not os.path.isfile('{}/config/indices.json'.format(famdir)): raise IOError('Missing indices.json: Deuterocol2 directory structure not complete for {}'.format(famdir))

			if not os.path.isdir('{}/../cut_pdbs'.format(famdir)):
				os.mkdir('{}/../cut_pdbs'.format(famdir))
			if VERBOSITY: info('Cutting PDBs for {}'.format(famdir))
			if not self.skip_cut: self.cut_pdbs(famdir)

		for famdir in self.famdirs:
			lockfn = '{}/.lockfile'.format(famdir)
			skip = False
			try: 
				if os.path.isfile(lockfn): 
					with open(lockfn) as f:
						warn('Found lockfile in {} dated {}, skipping'.format(famdir, f.read().strip()))
						skip = True
						continue
				else:
					with open(lockfn, 'w') as f: f.write(time.strftime('%Y-%m-%d %H:%M:%S'))
				self.main(famdir)
			finally: 
				if os.path.isfile(lockfn) and not skip: os.remove(lockfn)
		info('Finished all assigned alignments')


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--skip-cut', action='store_true')
	parser.add_argument('d2dir', default='deuterocol2', help='Deuterocol2 directory (contains config/, pdbs/, and superpositions/)')

	args = parser.parse_args()


	x = TMalign(d2dir=args.d2dir, skip_cut=args.skip_cut)
	x.run()
