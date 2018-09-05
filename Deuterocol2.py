#!/usr/bin/env python2

from __future__ import print_function
import subprocess, os, re, sys
import json
import argparse

import Deuterocol1

class Paragraph(object):
	def __init__(self, tmdatadir, d1dir, outdir, pdblist, force=False):
		self.tmdatadir = tmdatadir
		self.d1dir = d1dir
		self.outdir = outdir
		self.force = force

		self.fams1, self.fams2 = pdblist
		self.fam1 = sorted(self.fams1)[0]
		self.fam2 = sorted(self.fams2)[0]

		self.tcmap = {}

	@staticmethod
	def load_d2(d2obj, pdblist, prefix=None):
		selfobj = Paragraph(tmdatadir=d2obj.tmdatadir, d1dir=d2obj.d1dir, outdir=d2obj.outdir, pdblist=pdblist, force=d2obj.force)

		#selfobj.fams1, selfobj.fams2 = pdblist
		
		selfobj.outdir = '{}/{}_vs_{}'.format(selfobj.outdir, selfobj.fam1, selfobj.fam2)

		return selfobj

	def initialize_dir(self):
		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)

		for subdir in ('config', 'html', 'pdbs', 'sequences', 'superpositions'):
			if not os.path.isdir('{}/{}'.format(self.outdir, subdir)): 
				os.mkdir('{}/{}'.format(self.outdir, subdir))
		with open('{}/config/command_line'.format(self.outdir), 'w') as f: 
			f.write(str(sys.argv))

		with open('{}/config/align_me.json'.format(self.outdir), 'w') as f: 
			f.write(json.dumps(self.fams1) + '\n')
			f.write(json.dumps(self.fams2) + '\n')

		self.tcmap = self.get_tcmap()
		with open('{}/config/tcmap.json'.format(self.outdir), 'w') as f: f.write(json.dumps(self.tcmap))

	def get_tcmap(self):
		tcmap = {}
		pdbs = set()
		tcids = set()
		#for l in pdblist:
		#	for fam in l:
		#		tcids.add(fam)
		#		for pdb in l[fam]: pdbs.add(pdb)
		pdblist = {self.fam1:self.fams1[self.fam1], self.fam2:self.fams2[self.fam2]}
		for fam in pdblist: 
			tcids.add(fam)
			for pdb in pdblist[fam]: pdbs.add(pdb)
		pdbs = sorted(pdbs)
		tcids = sorted(tcids)
		with open('{}/tcmap.tsv'.format(self.tmdatadir)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.startswith('#'): continue
				else:
					for fam in tcids:
						#check if possibly relevant
						if fam in l and Deuterocol1.TCID.parse_str(l.split()[0]) in Deuterocol1.TCID.parse_str(fam):
							#see which PDB it is
							for pdb in pdbs:
								if pdb in l: 
									tcmap[pdb] = l.split()[0]
		return tcmap

class Deuterocol2(object):
	def __init__(self, tmdatadir, d1dir, outdir, force=False):
		self.tmdatadir = tmdatadir
		self.d1dir = d1dir
		self.outdir = outdir

		self.force = force

	def run(self, famlist1, famlist2):
		pdblist = []
		for famlist in (famlist1, famlist2):
			pdblist.append({})
			with open('{}/pdblist'.format(self.d1dir)) as f:
				for l in f:
					sl = l.split('\t')
					ttcid = Deuterocol1.TCID.parse_str(sl[0])
					for fam in famlist:
						qtcid = Deuterocol1.TCID.parse_str(fam)
						if ttcid in qtcid: 
							try: pdblist[-1][fam] += sl[1].strip().split(',')
							except KeyError: pdblist[-1][fam] = sl[1].strip().split(',')
		
		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)

		for fam1 in sorted(pdblist[0]):
			for fam2 in sorted(pdblist[1]):
				fam2pdb = [{fam1:pdblist[0][fam1]}, {fam2:pdblist[1][fam2]}]
				x = Paragraph.load_d2(d2obj=self, pdblist=fam2pdb)
				x.initialize_dir()
				print(pdblist)

				exit(1)


						

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--fams1', nargs='+', help='First list of families')
	parser.add_argument('--fams2', nargs='+', help='Second list of families')

	parser.add_argument('--d1dir', default='deuterocol1', help='Directory containing Deuterocol1 output')

	parser.add_argument('--tmdatadir', default='tmdata', help='Directory containing TM-prediction data and PDBs')
	parser.add_argument('--outdir', default='deuterocol2', help='Directory intended to contain Deuterocol2 output')

	args = parser.parse_args()

	if not (args.fams1 and args.fams2): 
		parser.print_usage()
		exit(1)

	deut = Deuterocol2(tmdatadir=args.tmdatadir, d1dir=args.d1dir, outdir=args.outdir)
	deut.run(args.fams1, args.fams2)
