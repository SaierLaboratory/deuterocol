#!/usr/bin/env python2
from __future__ import print_function, division

import sys
import os, re
import urllib
import argparse
import xml.etree.ElementTree as ET

COLORS = {'1':'\033[31m', '2':'\033[34m', 'B':'\033[33m', 'H':'\033[33m', 'F':'\033[32m', 'U':'\033[37m', 'I':'\033[35m', 'C':'\033[33m'}
RESET = '\033[0m'
TRANSMEMBRANE = ['H', 'C', 'B']

def warn(*msgs):
	for x in msgs: print('[WARNING]:', x, file=sys.stderr)

class PDBTM:
	def __init__(self, filename):
		#self.tree = ET.parse(filename)
		#self.root = self.tree.getroot()
		def strsum(l):
			s = ''
			for x in l: s += x.rstrip() + '\n'
			return s
		f = open(filename)
		s = []
		for l in f: s.append(l)
		#s = strsum(s[1:-1]).strip()
		s = strsum(s).strip()

		self.root = ET.fromstring(s)
		print(root)

def get_database(prefix='.'):
	if not prefix.endswith('/'): prefix += '/'
	print('Fetching database...', file=sys.stderr)
	db = urllib.urlopen('http://pdbtm.enzim.hu/data/pdbtmall')
	print('Saving database...', file=sys.stderr)
	f = open('%s/pdbtmall' % prefix, 'w')
	#for l in db: f.write(l)
	f.write(db.read())
	db.close()
	f.close()

def build_database(fn, prefix):
	print('Unpacking database...', file=sys.stderr)
	f = open(fn)
	db = f.read()
	f.close()
	firstline = 1
	header = ''
	entries = []
	pdbids = []
	for l in db.split('\n'):
		if firstline: 
			header += l
			firstline -= 1
			continue
		if 'PDBTM>' in l: continue
		if l.startswith('<?'): continue
		if l.startswith('<pdbtm'):
			a = l.find('ID=') + 4
			b = a + 4
			pdbids.append(l[a:b])
			entries.append(header)
		entries[-1] += l + '\n'
	if not prefix.endswith('/'): prefix += '/'
	if not os.path.isdir(prefix): os.mkdir(prefix)
	for entry in zip(pdbids, entries):
		f = open(prefix + entry[0] + '.xml', 'w')
		f.write(entry[1])
		f.close()

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

	def load_xml(self, element):
		self.label = element.attrib['CHAINID']
		for section in element: 
			if section.tag.endswith('SEQ'): 
				self.seq = re.sub('[^ACDEFGHIKLMNPQRSTVWY]', 'X', re.sub('\s+', '', section.text))
				if re.match('X+$', self.seq): self.seq = ''
			elif section.tag.endswith('REGION'):
				self.regions.append(Region(parent=self))
				self.regions[-1].load_xml(section)
		if self.label in self.parent.fakechains: self.fake = True

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
