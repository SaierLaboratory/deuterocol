#!/usr/bin/env python2

from __future__ import print_function
import subprocess, os, re, sys
import argparse

import Deuterocol1

class Deuterocol2(object):
	def __init__(self, tmdatadir, d1dir, outdir):
		self.tmdatadir = tmdatadir
		self.d1dir = d1dir
		self.outdir = outdir


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

		print(pdblist)

		if not os.path.isdir(args.outdir): os.mkdir(args.outdir)

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
