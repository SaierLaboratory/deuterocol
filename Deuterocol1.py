#!/usr/bin/env python2

from __future__ import print_function

import argparse, os, re, tempfile, subprocess
import urllib
import sys

import Bio.Blast.NCBIXML

#FIXME: expose this to args
CORES = 8

VERBOSITY = 1

def info(*things): print('[INFO]:', *things, file=sys.stderr)

def probably_nucleic(s):
	if s.count('T') + s.count('C') + s.count('U') + s.count('A') + s.count('G') >= (0.9 * len(s)): return True
	else: return False

def pdb_fetch_seq(pdblist):
	post_data = 'structureIdList='
	for pdbid in sorted(pdblist): post_data += '{},'.format(pdbid)
	post_data = post_data[:-1] + '&compressionType=uncompressed'

	url = urllib.urlopen('https://www.rcsb.org/pdb/download/viewFastaFiles.do', data=post_data)
	#with open('post_data', 'w') as f: f.write(post_data)
	rawseq = url.read()
	url.close()

	fastas = []
	for l in rawseq.split('\n'):
		if not l.strip(): continue
		elif l.startswith('>'): fastas.append(Fasta(header=l.strip()))
		else: fastas[-1].seq += l.strip()

	#for fasta in fastas: print(fasta)
	return fastas

#TODO: offload to common stuff
class TCID(object):
	def __init__(self):
		self.tcid = ['', '', '', '', '', '']
		self.tc_class, self.tc_subclass, self.tc_family, self.tc_subfamily, self.tc_transporter, self.tc_id = self.tcid

	@staticmethod
	def parse_str(s):
		ss = re.split('\.|-', s)
		tcid = TCID()
		for i, x in enumerate(ss):
			tcid.tcid[i] = x
		return tcid

	def __len__(self):
		n = 0
		for x in self.tcid: 
			if x: n += 1
			else: break
		return n

	def __iter__(self): return iter(self.tcid)

	def __hash__(self): return hash(tuple(self.tcid))

	def __contains__(self, other):
		if type(other) is str: return TCID.parse_str(other) in self
		else:
			if len(other) < len(self): 
				return False
			else:
				for x1, x2 in zip(other, self):
					if (x1 and x2):
						if x1 != x2: return False
					elif x2 and not x1: return False
					elif x1 and not x2: return True
					else: return True
				return True

	def __str__(self):
		s = ''
		for i, x in enumerate(self.tcid):
			if not x: break
			elif i == 0: s += self.tcid[0]
			elif i < 5: s += '.{}'.format(x)
			elif i == 5: s += '-{}'.format(x)
		return s

	def __repr__(self):
		return 'TCID({})'.format(self)

#TODO: offload to common stuff
class Fasta(object):
	def __init__(self, header='untitled', seq=''):
		self.header = header
		self.seq = seq

	def copy(self): return Fasta(header=self.header, seq=self.seq)
	def __len__(self): return len(re.sub('\s', '', self.seq.strip()))
	def __str__(self): return '{}\n{}'.format(self.header, self.seq)
	def __iter__(self): return iter(self.seq)
	def __hash__(self): return hash(self.seq)

#TODO: Unify container classes
class Subunit(object):
	def __init__(self, pdbid, letter, spans=[]):
		self.pdbid = pdbid
		self.letter = letter
		self.spans = spans

	def __str__(self):
		out = '{}_{}\t'.format(self.pdbid, self.letter)
		for span in self.spans:
			out += '{}-{},'.format(span[0], span[1])
		return out[:-1]

	@staticmethod
	def parse_str(s):
		ls = s.strip().split('\t')
		pdbid, letter = ls[0][:4], ls[0][5:]
		spans = []
		for spanstr in ls[1].split(','):
			if spanstr.startswith('-'): 
				start = 0
				end = spanstr[spanstr[1:].find('-'):]
			else:
				start, end = spanstr.split('-')
			spans.append((int(start), int(end)))
		return Subunit(pdbid, letter, spans)
			

class TMData(object):
	def __init__(self):
		self.subunits = {}

	def get_all_subunits(self):
		out = []
		for mode in self.subunits:
			[out.append('{}_{}'.format(subunit.pdbid, subunit.letter)) for subunit in self.subunits[mode]]
		return out

	def get_distinct_subunits(self):
		out = set()
		for mode in self.subunits:
			[out.add('{}_{}'.format(subunit.pdbid, subunit.letter)) for subunit in self.subunits[mode]]
		return sorted(out)

	def load_from_dir(self, indir):
		fns = []

		if VERBOSITY: info('Loading TM-assignment data...')
		for subdir in os.listdir(indir):
			if os.path.isdir('{}/{}'.format(indir, subdir)):
				self.subunits[subdir] = []
				if not os.path.isfile('{}/{}/ASSIGNMENTS.TSV'.format(indir, subdir)): continue
				fn = '{}/{}/ASSIGNMENTS.TSV'.format(indir, subdir)
				with open(fn) as f:
					for l in f:
						if not l.strip(): continue
						elif l.lstrip().startswith('#'): continue
						else: self.subunits[subdir].append(Subunit.parse_str(l))		

#TODO: let's try using sqlite instead of scattered TSVs
class Deuterocol1(object):
	def __init__(self, tmdatadir, outdir):
		self.tmdatadir = tmdatadir
		self.outdir = outdir
		self.tmdata = TMData()

	def get_pdb_sequences(self):
		subunitlist = self.tmdata.get_distinct_subunits()
		subunitlist = [x[:4].upper() + x[4:] for x in subunitlist]
		#print(subunitlist)

		tf = tempfile.NamedTemporaryFile()
		s = ''
		for pdbc in subunitlist: tf.write('{}\n'.format(pdbc))
		tf.flush()

		if VERBOSITY: info('Querying pdbaa for PDB sequences...')
		p = subprocess.Popen(['blastdbcmd', '-db', 'pdbaa', '-entry_batch', tf.name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

		out, err = p.communicate()

		#for pdbc in subunitlist
		direct = {}
		record = 0
		current = ''
		previous = ''
		print(len(subunitlist))
		sulist = subunitlist[:]

		fastas = []
		for l in out.split('\n'):
			if not l.strip(): continue
			elif l.startswith('>'):
				fastas.append(Fasta(header=l.strip()))
			else: fastas[-1].seq += l.strip()

		notfound = []
		for l in err.split('\n'):
			if not l.strip(): continue
			elif 'Entry not found' in l: notfound.append(l.split()[-1])

		for current in subunitlist:
			if current not in notfound:
				direct[current] = fastas.pop(0)
				newheader = direct[current].header
				newheader = newheader[newheader.find(current)-1:]
				newheader = newheader[:newheader.find('>', 1)]
				direct[current].header = newheader
			elif current[:4] == previous[:4] and previous in direct:
				direct[current] = direct[previous].copy()
				newheader = direct[current].header.replace(previous, current)
				direct[current].header = newheader
			else: pass
			previous = current

		#print(len(direct), len(notfound))

		#for sequences not found in PDBAA, try downloading from PDB
		fetchme = set()
		if VERBOSITY: info('Downloading sequences not found in pdbaa..')
		for current in subunitlist:
			if current not in direct: fetchme.add(current[:4])
		fastas = pdb_fetch_seq(fetchme)

		foundit = 0
		for fasta in fastas:
			fh = fasta.header[1:fasta.header.find('|')].replace(':', '_')

			for current in subunitlist:
				if current in direct: continue
				#elif current == fh:
				#	direct[current] = fasta.copy()
				#	direct[current].header = '>{}'.format(current)
				#desperation strikes
				elif current[:4] == fh[:4]:
					direct[current] = fasta.copy()
					direct[current].header = '>{}'.format(current)
				previous = current

		#print(len(direct))
		#at this point, 300/15208 (2%) of the targets are still missing somehow
		#however, most of them seem to be obsolete or theoretical models

		if VERBOSITY: info('Saving sequences...')
		if not os.path.isdir('{}/sequences'.format(self.tmdatadir)): 
			os.mkdir('{}/sequences'.format(self.tmdatadir))
		for subunit in direct:
			with open('{}/sequences/{}.fa'.format(self.tmdatadir, subunit), 'w') as f:
				x = direct[subunit].copy()
				x.header = '>{}|PDBID|CHAIN|SEQUENCE'.format(subunit.replace('_', ':'))
				f.write(str(x))

	def blast_against_tcdb(self):
		if VERBOSITY: info('BLASTing against TCDB...')
		megafasta = ''
		for fn in os.listdir('{}/sequences'.format(self.tmdatadir)):
			if fn.startswith('.') or not fn.endswith('fa'): pass
			else: 
				with open('{}/sequences/{}'.format(self.tmdatadir, fn)) as f:
					s = f.read()
					if probably_nucleic(s[s.find('\n'):]): continue
					megafasta += s.strip() + '\n'
		p = subprocess.Popen(['blastp', '-db', 'tcdb', '-comp_based_stats', 'no', '-outfmt', '5', '-evalue', '1000', '-num_alignments', '1', '-out', '{}/tcblast.xml'.format(self.tmdatadir), '-num_threads', str(CORES)], stdin=subprocess.PIPE)
		p.communicate(input=megafasta)

	def get_tcmapping(self, expect=1e-12, identities=40):
		#note: an expect cut off of 1e-12 seems to work well, allowing more distant homologs and turning away fragments and noise
		if VERBOSITY: info('Getting best TCDB hits...')
		tcmap = {}
		tcids = set()
		with open('{}/tcblast.xml'.format(self.tmdatadir)) as f:
			root = Bio.Blast.NCBIXML.parse(f)


			for q in root: 
				qpdbc = q.query[:6].replace(':', '_')
				for aln in q.alignments:
					for hsp in aln.hsps:
						if hsp.identities < identities: continue
						elif hsp.expect >= expect: continue
						#elif hsp.expect <= 1e-13: continue
						ttcid = TCID.parse_str(aln.accession)
						#print(ttcid, qpdbc, ttcid.tcid, '{:016x}'.format(hash(ttcid)))
						tcids.add(ttcid)
						try: tcmap[str(ttcid)].append(qpdbc)
						except KeyError: tcmap[str(ttcid)] = [qpdbc]

		with open('{}/tcmap.tsv'.format(self.tmdatadir), 'w') as f:
			for tcid in sorted(tcmap):
				line = '{}\t'.format(tcid)
				for pdbc in tcmap[tcid]:
					line += '{},'.format(pdbc)
				f.write(line[:-1] + '\n')

	def get_pdbs(self, *tclist):
		if VERBOSITY: info('Getting PDB list...')
		tcmap = {}
		tcids = set()
		qtclist = [TCID.parse_str(tcstr2) for tcstr2 in tclist]
		pdbclist = set()
		pdbcdict = {}
		with open('{}/tcmap.tsv'.format(self.tmdatadir)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.lstrip().startswith('#'): continue
				sl = l.strip().split('\t')
				tcids.add(sl[0])
				tcmap[sl[0]] = sl[1].split(',')

				#for tcstr1 in sorted(tcmap):
				ttcid = TCID.parse_str(sl[0])
				for qtcid in qtclist:
					if ttcid in qtcid: 
						for pdbc in sl[1].split(','): 
							pdbclist.add(pdbc)
							try: pdbcdict[str(qtcid)].append(pdbc)
							except KeyError: pdbcdict[str(qtcid)] = [pdbc]

		with open('{}/pdblist'.format(self.outdir), 'w') as f:
			for qtc in sorted(pdbcdict): 
				s = '{}\t'.format(qtc)
				for pdbid in pdbcdict[qtc]: s += '{},'.format(pdbid)
				f.write('{}\n'.format(s[:-1]))
		exit(1)

	def fetch_pdbs(self, force=False):
		pdbidlist = set()
		with open('{}/pdblist'.format(self.outdir)) as f:
			for l in f: pdbidlist.add(l[:4])

		tf = tempfile.NamedTemporaryFile()
		for pdbid in sorted(pdbidlist):
			if not force:
				if os.path.isfile('{}/{}.pdb.gz'): continue
				if os.path.isfile('{}/{}.pdb'): continue
			tf.write('https://files.rcsb.org/download/{}.pdb.gz\n'.format(pdbid))
		tf.flush()

		if not os.path.isdir('{}/pdbs'.format(self.outdir)): os.mkdir('{}/pdbs'.format(self.outdir))

		if VERBOSITY: info('Downloading PDBs...')
		cmd = ['wget', '-P', '{}/pdbs'.format(self.outdir), '-i', tf.name, '--no-check-certificate']
		if not force: cmd.append('-nc')
		subprocess.call(cmd)

		for pdbid in sorted(pdbidlist):
			if os.path.isfile('{}/pdbs/{}.pdb.gz'.format(self.outdir, pdbid)):
				subprocess.call(['gunzip', '{}/pdbs/{}.pdb.gz'.format(self.outdir, pdbid)])


	def run(self, *tclist, **kwargs):
		force = kwargs['force'] if 'force' in kwargs else False

		if not os.path.isdir(self.tmdatadir): os.mkdir(self.tmdatadir)
		self.tmdata.load_from_dir(self.tmdatadir)
		if force or not (os.listdir('{}/sequences'.format(self.tmdatadir))): self.get_pdb_sequences()
		if force or not (os.path.isfile('{}/tcblast.xml'.format(self.tmdatadir))): self.blast_against_tcdb()
		if force or not (os.path.isfile('{}/tcmap.tsv'.format(self.tmdatadir))): self.get_tcmapping()

		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)
		#if force or not (os.path.isfile('{}/pdblist'.format(self.outdir))): self.get_pdbs(*tclist)
		self.get_pdbs(*tclist)
		self.fetch_pdbs(force=force)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Cross-checks OPM and PDBTM against TCDB')

	parser.add_argument('tmdatadir', nargs='?', default='tmdata', help='Directory containing all TMS definitions, alpha and beta. Relevant files must be named ASSIGNMENTS.TSV.')
	parser.add_argument('-f', action='store_true', help='Enables clobbering')
	parser.add_argument('--fams', nargs='+', help='List of families to pick up')
	parser.add_argument('-o', '--outdir', default='deuterocol1', help='Directory to send output to')

	args = parser.parse_args()

	deut = Deuterocol1(tmdatadir=args.tmdatadir, outdir=args.outdir)
	deut.run(*args.fams, force=args.f)
