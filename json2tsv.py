#!/usr/bin/env python
from __future__ import print_function
import json, argparse, sys

def pretty(obj, k):
	if 'cov' in k: return '{:0.2%}'.format(obj[k])
	elif 'present' in k: return '{}-{}'.format(obj[k][0][0], obj[k][0][1])
	else: return str(obj[k])

def main(infile, outfile):
	fout = open(outfile, 'w')
	with open(infile) as fin:
		header = ''
		blacklist = set(['distances', 'qaligned', 'saligned', 'queryfn', 'subjectfn'])
		#length	matrix	qhel	qpresent	quality	query	rmsd	shel	spresensubject
		whitelist = ['rmsd', 'quality', 'length', 'tmcov', 'fullcov', 'qtmlen', 'qlen', 'stmlen', 'slen', 'qpresent', 'spresent']
		for l in fin:
			name, rawjson = l.split('\t')
			query, qchain, qhel, vs, subject, schain, shel = name.split('_')
			obj = json.loads(rawjson)

			if not header:
				header += '#name\t'
				header = '#name\tquery\tqchain\tqtms\tsubject\tschain\tstms\t'
				for k in whitelist:
					if k in blacklist: continue
					header += k + '\t'
				fout.write(header.strip() + '\n')
			s = ''
			s += '{}\t'.format(name)
			s += '{}\t{}\t{}\t'.format(query, qchain, qhel)
			s += '{}\t{}\t{}\t'.format(subject, schain, shel)
			for k in whitelist:
				if k in blacklist: continue
				try: s += pretty(obj, k) + '\t'
				except KeyError: s += '\t'
			fout.write(s.strip() + '\n')
			
			

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', nargs='?', default='/dev/stdin', help='input file (default:stdin)')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout', help='output file (default:stdout)')

	args = parser.parse_args()
	if args.infile == '/dev/stdin': 
		print('Now reading input from stdin...', file=sys.stderr)

	main(args.infile, args.outfile)
