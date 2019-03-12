#!/usr/bin/env python

from __future__ import print_function, division

import sys, os, argparse, json
import extendomatic

def main(infile, outfile):
	g = open(outfile, 'w')
	g.write('#RMSD\tCoverage\tTM-score\n')
	with open(infile) as f:
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue
			name, jstr = l.split('\t')
			data = extendomatic.unpack_obj(json.loads(jstr))
			g.write('{rmsd}\t{mincov}\t{quality}\n'.format(**data))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout')

	args = parser.parse_args()

	main(args.infile, args.outfile)
