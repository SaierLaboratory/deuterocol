#!/usr/bin/env python2

from __future__ import print_function
import subprocess, os, re, sys
import json, shutil
import argparse

import Deuterocol1

class Paragraph(object):
	def __init__(self, tmdatadir, d1dir, outdir, pdblist, force=False, bundle=4, loopless=False):
		self.tmdatadir = tmdatadir
		self.d1dir = d1dir
		self.outdir = outdir
		self.force = force
		self.bundle = bundle
		self.loopless = False

		self.fams1, self.fams2 = pdblist
		self.fam1 = sorted(self.fams1)[0]
		self.fam2 = sorted(self.fams2)[0]

		self.tcmap = {}

	@staticmethod
	def load_d2(d2obj, pdblist, prefix=None):
		selfobj = Paragraph(tmdatadir=d2obj.tmdatadir, d1dir=d2obj.d1dir, outdir=d2obj.outdir, pdblist=pdblist, force=d2obj.force, bundle=d2obj.bundle, loopless=d2obj.loopless)

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

		indices = {}
		with open('{}/indices.json'.format(self.d1dir)) as f: 
			indobj = json.loads(f.read())
			#for pdbid in sorted(indobj): indices[pdbid] = Deuterocol1.Spans.parse_json(indobj[pdbid])
			for pdbid in sorted(indobj): 
				indices[pdbid] = Deuterocol1.Spans()
				for span in indobj[pdbid]: indices[pdbid].add(span)

		if self.loopless: sourcedir = '{}/pdbs_loopless'.format(self.d1dir)
		else: sourcedir = '{}/pdbs'.format(self.d1dir)

		for fam in self.fams1:
			for pdb in self.fams1[fam]:
				shutil.copy('{}/{}.pdb'.format(sourcedir, pdb[:4]), '{}/pdbs/{}.pdb'.format(self.outdir, pdb[:4]))
		for fam in self.fams2:
			for pdb in self.fams2[fam]:
				shutil.copy('{}/{}.pdb'.format(sourcedir, pdb[:4]), '{}/pdbs/{}.pdb'.format(self.outdir, pdb[:4]))

		commands = []
		for fam1 in self.fams1:
			for fam2 in self.fams2:
				for pdb1 in self.fams1[fam1]:
					for pdb2 in self.fams2[fam2]:
						#XXX not sure if this works properly for short proteins
						if len(indices[pdb1]) <= self.bundle: endtms1 = 1
						else: endtms1 = len(indices[pdb1]) - self.bundle  + 1
						for bundle1 in range(0, endtms1):
							if len(indices[pdb2]) <= self.bundle: endtms2 = 1
							else: endtms2 = len(indices[pdb2]) - self.bundle + 1
							for bundle2 in range(0, endtms2):
								#print(fam1, pdb1, '{}-{}'.format(bundle1+1, bundle1+self.bundle-1+1), pdb2, '{}-{}'.format(bundle2+1, bundle2+self.bundle-1+1), fam2)
								qstart = indices[pdb1][bundle1].start
								if bundle1 + self.bundle - 1 >= len(indices[pdb1]): qend = indices[pdb1][-1].end
								else: qend = indices[pdb1][bundle1 + self.bundle - 1].end

								sstart = indices[pdb2][bundle2].start
								if bundle2 + self.bundle - 1 >= len(indices[pdb2]): send = indices[pdb2][-1].end
								else: send = indices[pdb2][bundle2 + self.bundle - 1].end
								commands.append({'name':'{}_h{}-{}_vs_{}_h{}-{}'.format( \
										pdb1, \
										bundle1+1, bundle1+self.bundle-1+1, \
										pdb2, \
										bundle2+1, bundle2+self.bundle-1+1, \
									), \
									'query': pdb1[:4], \
									'subject': pdb2[:4], \
									#'qhelices': list(range(bundle1+1, bundle1+1+self.bundle)), \
									#'shelices': list(range(bundle2+1, bundle2+1+self.bundle)), \
									'qhelices': [bundle1+1, bundle1+self.bundle], \
									'shelices': [bundle2+1, bundle2+self.bundle], \
									#'qfam': fam1, \
									#'sfam': fam2, \
									'bundle': self.bundle, \
									#'loopless': int(self.loopless), \
									'qchain': pdb1[-1], \
									'schain': pdb2[-1], \
									'qspan': [qstart, qend], \
									'sspan': [sstart, send], \
								})
								#print(commands[-1]['name'], fam1, fam2, pdb1, commands[-1]['qchain'], pdb2, commands[-1]['schain'])
		with open('{}/config/agenda.json'.format(self.outdir), 'w') as f: 
			for l in commands: f.write(json.dumps(l) + '\n')

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
	def __init__(self, tmdatadir, d1dir, outdir, force=False, bundle=4, loopless=False, allow_internal=False):
		self.tmdatadir = tmdatadir
		self.d1dir = d1dir
		self.outdir = outdir
		self.loopless = loopless
		self.allow_internal = allow_internal

		self.force = force

		self.bundle = bundle

	def run(self, famlist1, famlist2):
		if famlist2 == ['auto']:
			famlist2 = []
			with open('{}/pdblist.json'.format(self.d1dir)) as f:
				obj = json.loads(f.read())
				for tcid in obj:
					found = True
					for fam in famlist1:
						if Deuterocol1.parse_str(tcid) in Deuterocol1.parse_str(fam): 
							found = False
							break
						if not found: break
					if found: famlist2.append(tcid)
		pdblist = []
		for famlist in (famlist1, famlist2):
			pdblist.append({})
			with open('{}/pdblist.json'.format(self.d1dir)) as f:
				pdblistobj = json.loads(f.read())
				for tcid in pdblistobj: 
					ttcid = Deuterocol1.TCID.parse_str(tcid)
					for fam in famlist:
						qtcid = Deuterocol1.TCID.parse_str(fam)
						if ttcid in qtcid:
							try: pdblist[-1][fam] += pdblistobj[tcid]
							except KeyError: pdblist[-1][fam] = pdblistobj[tcid]
				#for l in f:
				#	sl = l.split('\t')
				#	ttcid = Deuterocol1.TCID.parse_str(sl[0])
				#	for fam in famlist:
				#		qtcid = Deuterocol1.TCID.parse_str(fam)
				#		if ttcid in qtcid: 
				#			try: pdblist[-1][fam] += sl[1].strip().split(',')
				#			except KeyError: pdblist[-1][fam] = sl[1].strip().split(',')
		
		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)

		done = []
		for fam1 in sorted(pdblist[0]):
			for fam2 in sorted(pdblist[1]):
				if (fam1, fam2) in done: continue
				else: 
					done.append((fam1,fam2))
					done.append((fam2,fam1))
				if not self.allow_internal and fam1 == fam2: continue
				fam2pdb = [{fam1:pdblist[0][fam1]}, {fam2:pdblist[1][fam2]}]
				x = Paragraph.load_d2(d2obj=self, pdblist=fam2pdb)
				x.initialize_dir()


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('-l', type=int, default=4, help='Bundle size')

	parser.add_argument('--fams1', nargs='+', help='First list of families')
	parser.add_argument('--fams2', nargs='+', help='Second list of families')

	parser.add_argument('--d1dir', default='deuterocol1', help='Directory containing Deuterocol1 output')

	parser.add_argument('--loopless', action='store_true', help='Trim off loops')

	parser.add_argument('--tmdatadir', default='tmdata', help='Directory containing TM-prediction data and PDBs')
	parser.add_argument('--outdir', default='deuterocol2', help='Directory intended to contain Deuterocol2 output')
	parser.add_argument('--allow-internal', action='store_true', help='Allow self-vs-self comparisons')

	args = parser.parse_args()

	if not (args.fams1 and args.fams2): 
		parser.print_usage()
		exit(1)

	deut = Deuterocol2(tmdatadir=args.tmdatadir, d1dir=args.d1dir, outdir=args.outdir, bundle=args.l, loopless=args.loopless, allow_internal=args.allow_internal)
	deut.run(args.fams1, args.fams2)
