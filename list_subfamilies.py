#!/usr/bin/env python

from __future__ import print_function, division, generators

import argparse, Deuterocol1

def main(tcids, depth=4, tmdatadir='tmdata'):
	entries = set()
	with open('{}/tcmap.tsv'.format(tmdatadir)) as f:
		for l in f:
			if tcids:
				for tcid in tcids:
					if Deuterocol1.TCID.parse_str(l.split('\t')[0]) in Deuterocol1.TCID.parse_str(tcid):
						entries.add(str(Deuterocol1.TCID.parse_str(l.split('\t')[0])[:depth]))
			else:
				entries.add(str(Deuterocol1.TCID.parse_str(l.split('\t')[0])[:depth]))
	for e in sorted(entries): print(e)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('--tmdatadir', default='tmdata', help='where the tmdata directory is')
	parser.add_argument('-d', type=int, default=4, help='level of subset (1: class, 2: subclass, 3: family/superfamily, 4: subfamily/family, 5: system, 6: component, 7: PDB')
	parser.add_argument('tcid', nargs='*', default=[], help='list of TCIDs to get subsets for')

	args = parser.parse_args()
	if (args.d < 1) or (args.d > 7): 
		print('[ERROR]: Invalid depth')
		parser.print_help()
		exit(1)

	main(args.tcid, depth=args.d, tmdatadir=args.tmdatadir)
