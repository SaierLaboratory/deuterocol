#!/usr/bin/env python2

from __future__ import print_function, division, generators

import argparse, json, os, subprocess, re, sys, time
import tempfile
import superpose
import shutil


VERBOSITY = 1
def info(*things):
	print('[INFO]:', *things, file=sys.stderr)
def warn(*things):
	print('[WARNING]:', *things, file=sys.stderr)
def error(*things):
	print('[ERROR]:', *things, file=sys.stderr)
	exit(1)

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
					else: done.append(l.split('\t')[0])
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
									aligned.append(tml.strip())
						sp.quality = max(tmscores)
						qseq, dseq, sseq = aligned
						dseq += ' ' * max(0, len(qseq)-len(dseq))
						qaligned, saligned, distances = [], [], []
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
								elif dist == ':':
									if qcur == '-': qaligned.append([None, 1])
									else: qaligned.append([qi, 1])
									if scur == '-': saligned.append([None, 1])
									else: saligned.append([si, 1])
									distances.append(0.0)
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
							elif dist == ':': 
								if qaligned[-1][0] is None: qaligned.append([qi, 1])
								else: qaligned[-1][1] += 1
								if saligned[-1][0] is None: saligned.append([si, 1])
								else: saligned[-1][1] += 1
								distances.append(0.)
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
							lastqr = qcur
							lastsr = scur
							lastdist = dist

						sp.qaligned = qaligned
						sp.saligned = saligned
						sp.distances = distances
							
						query, qchain, qhel, vs, subject, schain, shel = obj['name'].split('_')
						sp.qhel = [int(x) for x in qhel[1:].split('-')]
						sp.shel = [int(x) for x in shel[1:].split('-')]
						f.write('{}\t{}\n'.format(obj['name'], sp.dump_json()))
						n += 1
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

		for pdb in indices:
			if len(indices[pdb]) <= 2: continue
			if bundle >= len(indices[pdb]): 
				shutil.copy('{}/../pdbs/{}.pdb'.format(famdir, pdb[:4]), 
					'{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdb, 1, len(indices[pdb])))
			else:
				for start in range(0, len(indices[pdb])-bundle+1):
					end = start + bundle - 1
					infile = '{}/../pdbs/{}.pdb'.format(famdir, pdb[:4])
					outfile = '{}/../cut_pdbs/{}_h{}-{}.pdb'.format(famdir, pdb, start+1, end+1)
					if not os.path.isfile(infile): continue
					if os.path.isfile(outfile) and os.path.getsize(outfile):
						with open(outfile) as f:
							for l in f: pass
						if l.startswith('END'): continue
					cmd = 'lvresidue {}-{}\nlvchain {}\nwrite PDB'.format(
						indices[pdb][start][0],
						indices[pdb][end][1],
						pdb[-1])
					p = subprocess.Popen(['pdbcur', 
						'xyzin', infile,
						'xyzout', outfile],
						stdin=subprocess.PIPE,
						stdout=subprocess.PIPE,
						stderr=subprocess.PIPE)
					out, err = p.communicate(input=cmd)
					if VERBOSITY: print(out)


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
