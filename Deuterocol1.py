#!/usr/bin/env python2

from __future__ import print_function

import argparse, os, re, tempfile, subprocess, json
import urllib
import sys

import Bio.Blast.NCBIXML

#FIXME: expose this to args
CORES = 8

VERBOSITY = 1

def info(*things): print('[INFO]:', *things, file=sys.stderr)
def warn(*things): print('[WARNING]:', *things, file=sys.stderr)

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
class Interval(object):
	def __init__(self, a, b):
		self.start = a
		self.end = b
	def __contains__(self, something):
		return True if self.start <= something <= self.end else False

	def __iter__(self): return self.start, self.end
	def __repr__(self): return 'Interval({}, {})'.format(self.start, self.end)
	def intersects(self, other):
		if self.start <= other.start <= self.end: return True
		elif self.start <= other.end <= self.end: return True
		elif other.start <= self.start <= other.end: return True
		elif other.start <= self.end <= other.end: return True
		else: return False
	def union(self, other):
		return Interval(min(self.start, other.start), max(self.end, other.end))
	def __lt__(self, other):
		if self.start < other.start: return True
		elif self.start == other.start:
			if self.end < other.end: return True
			else: return False
		else: return False

class Spans(object):
	def __init__(self, spans=None):
		self.spans = [] if spans is None else spans

	@staticmethod
	def parse_str(s):
		spansobj = Spans()
		for ss in re.split('\s*,\s*', s.strip()):
			spansobj.add([int(i) for i in ss.split('-')])
		return spansobj

	@staticmethod
	def parse_json(s):
		obj = json.loads(s)
		intervals = []
		for x in obj: intervals.append(Interval(x[0], x[1]))
		return Spans(intervals)

	def dump_json(self):
		obj = [[s.start, s.end] for s in self]
		return json.dumps(obj)

	def add(self, span): self.spans.append(Interval(*span))

	def __len__(self): return len(self.spans)
	def __repr__(self): 
		s = 'Span(Interval[{}])'.format(len(self))
		return s
	def __iter__(self): return iter(self.spans)
	def __str__(self):
		s = 'Span(['
		#for i in self: s += str(i) + ', '
		for i in self: s += '{}-{}, '.format(i.start, i.end)
		s = s[:-2] + '])'
		return s
	def __getitem__(self, i): return self.spans[i]
	def __setitem__(self, i, x): self.spans[i] = x

	def extend(self, other, selfish=False):
		mergeme = []
		exc1 = set()
		exc2 = set()
		for i, s1 in enumerate(self):
			for j, s2 in enumerate(other):
				#if s1.intersects(s2): 
				#this is a quick, dirty workaround that will surely be replaced by geometric constraints
				if s1.intersects(s2) and max(s2.end - s1.start, s1.end - s2.start) <= 40:  
					mergeme.append((i,j))
					exc1.add(i)
					exc2.add(j)
				if selfish: exc2.add(j)
		spansobj = Spans()
		for s in mergeme: 
			newinterval = self[s[0]].union(other[s[1]])
			spansobj.add([newinterval.start, newinterval.end])
		for i, s in enumerate(self):
			if i in exc1: continue
			spansobj.add([self[i].start, self[i].end])
		for j, s in enumerate(other):
			if j in exc2: continue
			spansobj.add([other[j].start, other[j].end])
		spansobj.spans.sort()
		return spansobj

	def merge(self):
		unmerged = True
		while unmerged:
			unmerged = False
			for i in range(len(self)-1):
				#if self[i].intersects(self[i+1]):
				#this is a quick, dirty workaround that will surely be replaced by geometric constraints
				if self[i].intersects(self[i+1]) and max(self[i+1].end - self[i].start, self[i].end - self[i].start) <= 40:
					unmerged = True
					self[i] = self[i].union(self[i+1])
					self.spans.pop(i+1)
					break
		

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

		pdb2tcdb = {}
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
							pdb2tcdb[pdbc] = str(ttcid)
		with open('{}/tcmap.json'.format(self.outdir), 'w') as g: g.write(json.dumps(pdb2tcdb))

		with open('{}/pdblist.json'.format(self.outdir), 'w') as f: f.write(json.dumps(pdbcdict))

	def fetch_pdbs(self, force=False):
		pdbidlist = set()
		with open('{}/pdblist.json'.format(self.outdir)) as f:
			pdbcdict = json.loads(f.read())
			for fam in pdbcdict: 
				[pdbidlist.add(pdbid[:4]) for pdbid in pdbcdict[fam]]

		tf = tempfile.NamedTemporaryFile()
		dl = False
		for pdbid in sorted(pdbidlist):
			if not force:
				if os.path.isfile('{}/pdbs/{}.pdb.gz'.format(self.outdir, pdbid)): continue
				elif os.path.isfile('{}/pdbs/{}.pdb'.format(self.outdir, pdbid)): continue
			tf.write('https://files.rcsb.org/download/{}.pdb.gz\n'.format(pdbid))
			dl = True
		tf.flush()

		if not os.path.isdir('{}/pdbs'.format(self.outdir)): os.mkdir('{}/pdbs'.format(self.outdir))

		if dl:
			if VERBOSITY: info('Downloading PDBs...')
			cmd = ['wget', '-P', '{}/pdbs'.format(self.outdir), '-i', tf.name, '--no-check-certificate']
			if not force: cmd.append('-nc')
			else:
				for fn in os.listdir('{}/pdbs'.format(self.outdir)):
					if fn.startswith('.'): continue
					elif fn.endswith('.pdb'): os.remove('{}/pdbs/{}'.format(self.outdir, fn))
					elif fn.endswith('.pdb.gz'): os.remove('{}/pdbs/{}'.format(self.outdir, fn))
			if dl: subprocess.call(cmd)
		elif VERBOSITY: info('Skipped downloading PDBs')

		for pdbid in sorted(pdbidlist):
			if os.path.isfile('{}/pdbs/{}.pdb.gz'.format(self.outdir, pdbid)):
				cmd = ['gunzip']
				if force: cmd.append('-f')
				cmd.append('{}/pdbs/{}.pdb.gz'.format(self.outdir, pdbid))
				subprocess.call(cmd)

	def fetch_indices_method1(self):
		''' Method 1: OPM, PDBTM, and STRIDE

		$OPM $PDBTM union $STRIDE intersection
		'''
		pdbidlist = set()
		with open('{}/pdblist.json'.format(self.outdir)) as f:
			pdbcdict = json.loads(f.read())
			for fam in pdbcdict: 
				[pdbidlist.add(pdbid) for pdbid in pdbcdict[fam]]

		opmspans = self.fetch_spans(pdbidlist, 'opm')
		pdbtmspans = self.fetch_spans(pdbidlist, 'pdbtm')
		stridespans = self.gen_stride_spans(pdbidlist, 'stride')

		unionspans = {}
		tmspans = {}

		indextable = {}
		for pdbid in sorted(pdbidlist):
			unionspans[pdbid] = opmspans[pdbid].extend(pdbtmspans[pdbid])
			tmspans[pdbid] = unionspans[pdbid].extend(stridespans[pdbid], selfish=True)
			tmspans[pdbid].merge()

			indextable[pdbid] = [[s.start, s.end] for s in tmspans[pdbid]]
		with open('{}/indices.json'.format(self.outdir), 'w') as f:
			f.write(json.dumps(indextable))

			#print(pdbid)
			#print('\t', opmspans[pdbid])
			#print('\t', pdbtmspans[pdbid])
			#print('\t', unionspans[pdbid])
			#print('\t', stridespans[pdbid])
			#print('\t', tmspans[pdbid])
			#print('#'*80)

	def gen_stride_spans(self, pdbidlist, db='stride'):
		#TODO: a practical way to fold this back into a dbtool
		spans = {}
		altspans = {}
		for pdbid in sorted(pdbidlist): spans[pdbid] = Spans()
		pdbs = set([pdbid[:4] for pdbid in pdbidlist])
		for pdb in pdbs: altspans[pdb] = Spans()

		for pdb in pdbs:
			cmd = ['stride', '{}/pdbs/{}.pdb'.format(self.outdir, pdb)]
			out = subprocess.check_output(cmd)
			for l in out.split('\n'):
				if not l.strip(): continue
				#TODO: generalize to allow beta barrels too
				elif l.startswith('LOC  AlphaHelix'): 
					#print(l)
					start = int(l[21:27].strip())
					chain = l[28]
					end = int(l[39:45].strip())
					try: spans[pdb + '_' + chain].add([start, end])
					except KeyError: altspans[pdb].add([start, end])
		#this could probably be reworked
		for pdbid in spans:
			if len(spans[pdbid]) == 0:
				spans[pdbid] = altspans[pdbid[:4]]
		return spans
						

	def fetch_spans(self, pdbidlist, db='opm'):
		spans = {}
		altspans = {}
		for pdbid in sorted(pdbidlist): spans[pdbid] = Spans()
		for pdbid in sorted(pdbidlist): altspans[pdbid[:4]] = Spans()

		with open('{}/{}/ASSIGNMENTS.TSV'.format(self.tmdatadir, db)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.startswith('#'): continue
				sl = l.split()
				dbid = sl[0][:4].upper() + sl[0][4:]
				if dbid in spans:
					spans[dbid] = Spans.parse_str(sl[1])
				elif dbid[:4] in altspans:
					altspans[dbid[:4]] = Spans.parse_str(sl[1])
		for pdbid in spans:
			#if spans[pdbid] is None:
			if len(spans[pdbid]) == 0:
				if VERBOSITY: warn('PDB {} not listed on {}, falling back on {}'.format(pdbid, db, pdbid[:4]))
				spans[pdbid] = altspans[pdbid[:4]]
		return spans

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

		self.fetch_indices_method1()

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Cross-checks OPM and PDBTM against TCDB')

	parser.add_argument('tmdatadir', nargs='?', default='tmdata', help='Directory containing all TMS definitions, alpha and beta. Relevant files must be named ASSIGNMENTS.TSV.')
	parser.add_argument('-f', action='store_true', help='Enables clobbering')
	parser.add_argument('--fams', nargs='+', help='List of families to pick up')
	parser.add_argument('-o', '--outdir', default='deuterocol1', help='Directory to send output to')

	args = parser.parse_args()

	deut = Deuterocol1(tmdatadir=args.tmdatadir, outdir=args.outdir)
	deut.run(*args.fams, force=args.f)