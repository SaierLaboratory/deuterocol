#!/usr/bin/env python
from __future__ import print_function, division
import argparse
import numpy as np

def main(infile, outfile, count=None):
	with open(infile) as f:
		n = 0
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue
			n += 1
			if count is None: print(l.rstrip())
		if count is None: return 0

		f.seek(0)

		if count < n:
			grabme = [int(x) for x in np.arange(0, n, n/count)]
			for i, l in enumerate(f):
				if not grabme: break
				elif i == grabme[0]:
					grabme.pop(0)
					print(l.rstrip())
		else:
			for l in f:
				print(l.rstrip())


if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('-n', type=int, default=None)
	parser.add_argument('infile')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout')
	args = parser.parse_args()

	main(args.infile, args.outfile, count=args.n)
