#!/usr/bin/env python2

from __future__ import print_function

import argparse, os, re, tempfile, subprocess
import urllib
import sys

import Bio.Blast.NCBIXML

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
		self.fake = False

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

	@staticmethod
	def parse_pdbtm_xml(element, parent=None):
		letter = element.attrib['CHAINID']
		regions = []
		for section in element: 
			if section.tag.endswith('SEQ'): 
				seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', 'X', re.sub('\s+', '', section.text))
				if re.match('X+$', seq): seq = ''
			elif section.tag.endswith('REGION'):
				regions.append(Region(parent=self))
				regions[-1].load_xml(section)
		if parent is not None and letter in parent.fakechains: self.fake = True
		return Subunit(pdbid=parent.pdbid, letter=letter, spans=regions)

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

class PDB(object):
	def __init__(self):
		self.label = ''
		self.chains = []
		self.fakechains = []

	def load_xml(self, f): 
		tree = ET.parse(f)
		root = tree.getroot()
		self.label = root.attrib['ID']
		self.fakechains = []
		for section in root: 
			if section.tag.endswith('BIOMATRIX'):
				for matrix in section:
					if matrix.tag.endswith('DELETE'):
						self.fakechains.append(matrix.attrib['CHAINID'])
						continue
					for subsection in matrix:
						if subsection.tag.endswith('APPLY_TO_CHAIN'): 
							self.fakechains.append(subsection.attrib['NEW_CHAINID'])
			if section.tag.endswith('CHAIN'): 
				self.chains.append(Chain(parent=self))
				self.chains[-1].load_xml(section)
	def get_fasta(self):
		s = ''
		hdone = []
		for c in self:
			if c.fake: continue
			if hash(c) in hdone: continue
			if not len(c.seq): continue
			s += '{}\n'.format(c.get_fasta())
			hdone.append(hash(c))

			#100 chains isn't even possible
			if len(hdone) > 100: hdone.pop(0)
		return s.rstrip()
	def __iter__(self): return iter(self.chains)

	def __str__(self):
		s = ''
		for c in self: 
			if c.label in self.fakechains: continue
			elif c.count() == 0: continue
			else: s += '{}_{}\n'.format(self.label, c)
		return s.strip()

class Chain(object):
	def __init__(self, label='_', seq='', regions=[], parent=None):
		self.label = label
		self.seq = seq
		self.regions = []
		self.parent = parent
		self.fake = False

	def get_fasta(self):
		return '>pdb|{}|{}\n{}'.format(self.parent.label, self.label, self.seq)

	def __hash__(self): return hash(self.seq)

	def __iter__(self): return iter(self.regions)

	def count(self, topotypes='BCH'):
		s = 0
		for r in self.regions:
			if r.topo in topotypes: s += 1
		return s

	def __str__(self):
		s = ''
		#s += str(self.parent)
		#s += '\t'
		s += str(self.label)
		s += '\t'
		for r in self.regions:
			if r.topo in 'BCH': s += '{}-{},'.format(r.pdb_beg, r.pdb_end)
		return s[:-1]

	def save(self, outfile):
		pass

class Region(object):
	def __init__(self, label=-1, seq='', topo='X', parent=None):
		self.label = label
		self.seq = seq
		self.topo = topo
		self.parent = parent

		self.seq_beg, self.seq_end = None, None
		self.pdb_beg, self.pdb_end = None, None

	def load_xml(self, element): 
		self.seq_beg = int(element.attrib['seq_beg'])
		self.seq_end = int(element.attrib['seq_end'])
		self.pdb_beg = int(element.attrib['pdb_beg'])
		self.pdb_end = int(element.attrib['pdb_end'])
		self.topo = element.attrib['type']
		if self.parent is not None:
			self.seq = self.parent.seq[self.seq_beg-1:self.seq_end]

	def __contains__(self, i):
		#not sure if this should be seq_beg or pdb_beg
		if self.seq_beg-1 <= i <= self.seq_end: return True
		else: return True

	def __iter__(self): return iter(self.seq)

def extract_sequences(pdbtmdir):
	pdbs = []
	fout = open('{}/pdbtm.fasta'.format(pdbtmdir), 'w')
	for fn in os.listdir(pdbtmdir):
		if not re.match('[0-9][a-z]..\.xml', fn): continue
		with open('{}/{}'.format(pdbtmdir, fn)) as f:
			pdbs.append(PDB())
			pdbs[-1].load_xml(f)
			fout.write(pdbs[-1].get_fasta() + '\n')

def write_assignments(pdbtmdir):
	with open('{}/ASSIGNMENTS.TSV'.format(pdbtmdir), 'w') as g:
		for fn in os.listdir(pdbtmdir):
			if not re.match('[0-9][a-z]..\.xml', fn): continue
			with open('{}/{}'.format(pdbtmdir, fn)) as f:
				p = PDB()
				p.load_xml(f)
				g.write(str(p) + '\n')

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Manages PDBTM databases. Automatically fetches the PDBTM database if no options are specified. Run without any arguments, dbtool will retrieve the PDBTM database, store it in pdbtm, and unpack it.')

	parser.add_argument('-d', '--db', default='pdbtmall', help='name of concatenated database file {default:pdbtmall}')
	parser.add_argument('-b', '--build-db', action='store_true', help='(re)build database from an existing pdbtmsall file (available at http://pdbtm.enzim.hu/data/pdbtmall)')
	parser.add_argument('directory', nargs='?', default='pdbtm', help='directory to store database in')
	parser.add_argument('-f', '--force-refresh', action='store_true', help='force overwrite of existing database. Functionally equivalent to removing the old database and rerunning.')
	parser.add_argument('-o', '--outdir', default='tmdir', help='where to store TMS assignments')
	#parser.add_argument('-x', '--extract-sequences', action='store_true', help='extract sequences from PDBTM database too')
	#parser.add_argument('-n', metavar='bundle_size', type=int, help='size to cut bundles into')

	args = parser.parse_args()

	if args.build_db: build_database(args.db, args.directory)
	else: #db = PDBTM(args.db)
		if not os.path.isdir(args.directory): os.mkdir(args.directory)
		if args.force_refresh or not os.path.isfile('%s/%s' % (args.directory, args.db)): get_database(args.directory)
		build_database('%s/%s' % (args.directory, args.db), args.directory)

	#if args.extract_sequences:
	extract_sequences(args.directory)
	write_assignments(args.directory)
		

	#http://pdbtm.enzim.hu/data/pdbtmall
#!/usr/bin/env python
from __future__ import print_function, division

import urllib
import sys, os, re
import argparse
import json

#TODO: move common classes out
import pdbtm_dbtool

class OPM(object):
	def __init__(self):
		pass

def get_database(prefix='.'):
	print('Fetching database...', file=sys.stderr)
	db = urllib.urlopen('http://opm.phar.umich.edu/OPM-2018-07-26.json')
	print('Saving database...', file=sys.stderr)
	with open('{}/opmall.json'.format(prefix), 'w') as f:
		f.write(db.read())

def build_database(fn, prefix):
	print('Unpacking database...', file=sys.stderr)

	db = {'atlas':[], 'citation':[], 'class':[], 'classification':[], 'family':[], 'history':[], 'membrane':[], 'protein':[], 'relatedproteins':[], 'relatedproteinstemp':[], 'relatedproteinstempnew':[], 'species':[], 'subunits':[], 'superfamily':[], 'type':[]}
	with open(fn) as f:
		n = 0
		for l in f:
			if not l.startswith('['): continue
			db[sorted(db)[n]] = json.loads(l)
			n += 1

	proteins = {}
	for protein in db['protein']:
		pdbid = protein['pdbid']
		proteins[pdbid] = pdbtm_dbtool.PDB()
		proteins[pdbid].label = pdbid

	related_proteins = {}
	parent_proteins = {}
	for relprot in db['relatedproteins']:
		rel_pdbid = relprot['pdbid']
		par_pdbid = relprot['protein_pdbid']

		if par_pdbid in proteins and rel_pdbid in proteins: pass
		elif par_pdbid in proteins:
			parent_proteins[rel_pdbid] = par_pdbid
			try: related_proteins[par_pdbid].append(rel_pdbid)
			except KeyError: related_proteins[par_pdbid] = [rel_pdbid]
			#proteins[rel_pdbid] = pdbtm_dbtool.PDB()
			#proteins[rel_pdbid].label = rel_pdbid
		#necessary because 1ym6 is misannotated as being a relative of 1okc
		else:
			parent_proteins[par_pdbid] = rel_pdbid
			try: related_proteins[rel_pdbid].append(par_pdbid)
			except KeyError: related_proteins[rel_pdbid] = [par_pdbid]
			#proteins[par_pdbid] = pdbtm_dbtool.PDB()
			#proteins[par_pdbid].label = par_pdbid

	#print(related_proteins['1ym6'])
	#for p in related_proteins: print(p, related_proteins[p])
	#for p in proteins: print(p)
	#exit()

	print('Writing assignments to {}/ASSIGNMENTS.TSV...'.format(prefix), file=sys.stderr)
	writeme = []
	for subunit in db['subunits']:
		pdbid = subunit['pdbid']
		letter = subunit['letter']
		#print(subunit['pdbid'], subunit['letter'], subunit['pdbid'] in proteins)
		if pdbid not in proteins:
			#necessary because 2iqo_? is an orphaned subunit
			if pdbid in parent_proteins: continue
			else: proteins[pdbid] = pdbtm_dbtool.PDB()
		chain = pdbtm_dbtool.Chain(label=letter, parent=proteins[pdbid])

		proteins[pdbid].chains.append(chain)
		segments = subunit['segment']
		spans = re.split(r'\),?\s*[0-9][0-9]?\s*\(?', segments)
		segstr = ''
		trouble = '2he6'
		for i in spans:
			s = i
			if pdbid == trouble: print('start:', s)
			s = re.sub(r'\n[A-Za-z]+', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\r?\n', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'^\s*[0-9]+\(\s*', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\s+', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'[0-9]+\(', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\)[^\-]*', '', s)
			if pdbid == trouble: print(s)
			if s.strip(): segstr += s + ','
			if pdbid == trouble: print('final:', s)
		#print(segments)
		#print(segstr[:-1])
		#print(sorted(db))
		writeme.append('{}_{}\t{}\n'.format(pdbid, letter, segstr[:-1]))
		if pdbid == trouble: print(writeme[-1])
		if pdbid in related_proteins:
			for relpdbid in related_proteins[pdbid]:
				writeme.append('{}_{}\t{}\n'.format(relpdbid, letter, segstr[:-1]))

	last100 = []
	with open('{}/ASSIGNMENTS.TSV'.format(prefix), 'w') as f:
		for s in sorted(writeme):
			if not last100: pass
			elif hash(s) in last100: continue

			f.write(s)
			last100.append(hash(s))
			if len(last100) > 100: last100.pop(0)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Manages OPM databases. Automatically fetches the OPM database if no options are specified. Run without any arguments, opm_dbtool will retrieve the OPM database, store i in opm, and unpack it.')

	parser.add_argument('-d', '--db', default='opmall.json', help='name of concatenated database file (default:opmall.json)')
	parser.add_argument('-b', '--build-db', action='store_true', help='(re)build database from an existing opmall.json file (available at http://opm.phar.umich.edu/OPM-2018-07-26.json)')
	parser.add_argument('directory', nargs='?', default='opm', help='directory to store database in')
	parser.add_argument('-f', '--force-refresh', action='store_true', help='force overwrite of existing database. Functionally equivalent to removing the old database and rerunning.')
	parser.add_argument('-x', '--extract-sequences', action='store_true', help='extract sequences from OPM database too')

	args = parser.parse_args()

	if args.build_db: build_database(args.db, args.directory)
	else:
		if not os.path.isdir(args.directory): os.mkdir(args.directory)
		if args.force_refresh or not os.path.isfile('{}/{}'.format(args.directory, args.db)): get_database(args.directory)
		build_database('{}/{}'.format(args.directory, args.db), args.directory)

	#if args.extract_sequences:

