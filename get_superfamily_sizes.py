#!/usr/bin/env python

import json, os, argparse
import Deuterocol1

def main(tmdatadir, count=False):
	with open('{}/tcdb/superfamily.json'.format(tmdatadir)) as f:
		obj = json.loads(f.read())
	tcmap_chains = {}
	tcmap_strucs = {}
	with open('{}/tcmap.tsv'.format(tmdatadir)) as f:
		for l in f:
			sl = l.split('\t')
			tcmap_chains[sl[0]] = sl[1].split(',')
			tcmap_strucs[sl[0]] = set()
			[tcmap_strucs[sl[0]].add(pdbid[:4]) for pdbid in sl[1].split(',')]
			tcmap_strucs[sl[0]] = sorted(tcmap_strucs[sl[0]])

	tcids = tcmap_chains.keys()

	sfchains = {}
	sfstrucs = {}
	printme = []
	for superfam in obj:
		sfchains[superfam] = []
		sfstrucs[superfam] = []
		for fam in obj[superfam]:
			popme = []
			for i, tcid in enumerate(tcids):
				if Deuterocol1.TCID.parse_str(tcid) in Deuterocol1.TCID.parse_str(fam):
					popme.append(i)
					sfchains[superfam] += tcmap_chains[tcid]
					sfstrucs[superfam] += tcmap_strucs[tcid]
			for i in popme[::-1]: tcids.pop(i)
			#print(fam)
		printme.append((len(sfstrucs[superfam]), len(sfchains[superfam]), superfam))
	for strucs, chains, name in sorted(printme)[::-1]:
		print('{}: {} structures, {} chains'.format(name, strucs, chains))
	

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('tmdatadir', nargs='?', default='tmdata')
	parser.add_argument('-c', action='store_true', help='return counts')
	

	args = parser.parse_args()

	main(args.tmdatadir, count=args.c)
