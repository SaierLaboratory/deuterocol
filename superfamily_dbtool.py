#!/usr/bin/env python3

import re, json
import argparse

def parse_txtdump(f):
	superfamilies = {}
	lastsuperfam = ''

	mode = 2
	for l in f:
		if not l.strip(): continue
		print(l, mode)

		if mode == 2:
			if re.match('[0-9]+\.[A-Z]+\.[0-9]+ - ', l): 
				superfamilies[lastsuperfam].append(l.split()[0])
				continue
			elif re.search('[0-9]+\.[A-Z]+\.[0-9]+', l):
				superfamilies[lastsuperfam] += re.findall('[0-9]+\.[A-Z]+\.[0-9]+', l)
				continue
			elif len(l.split()) < 15 and l.count('.') < 5:
				mode = 0
				lastsuperfam = l.strip()
				superfamilies[lastsuperfam] = []
				continue

		elif mode == 1:
			if re.match('[0-9]+\.[A-Z]+\.[0-9]+ - ', l): 
				mode = 2
				superfamilies[lastsuperfam].append(l.split()[0])
				continue
			elif re.search('[0-9]+\.[A-Z]+\.[0-9]+', l):
				superfamilies[lastsuperfam] += re.findall('[0-9]+\.[A-Z]+\.[0-9]+', l)
				mode = 2
				continue
			elif 2 < len(l.split()) <= 15 and l.count('.') < 5:
				mode = 1
				lastsuperfam = l.strip()
				superfamilies[lastsuperfam] = []
				continue
			else: continue

		elif mode == 0:
			if re.match('[0-9]+\.[A-Z]+\.[0-9]+ - ', l): 
				mode = 2
				superfamilies[lastsuperfam].append(l.split()[0])
				continue
			elif re.search('[0-9]+\.[A-Z]+\.[0-9]+', l):
				superfamilies[lastsuperfam] += re.findall('[0-9]+\.[A-Z]+\.[0-9]+', l)
				mode = 2
				continue
			else: 
				mode = 1
				continue
	for superfam in superfamilies:
		superfamilies[superfam] = sorted(set(superfamilies[superfam]))
	return superfamilies
			
	

def main(infile, outfile, reverse=False):
	with open(infile) as f:
		obj = parse_txtdump(f)
		if reverse:
			robj = {}
			for superfam in obj:
				for fam in obj[superfam]:
					robj[fam] = superfam
			with open(outfile, 'w') as g: g.write(json.dumps(robj, indent=3, sort_keys=True))
		else:
			with open(outfile, 'w') as g: g.write(json.dumps(obj, indent=3, sort_keys=True))

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Preprocesses plaintext dumps of http://tcdb.org/superfamily.php into json')

	parser.add_argument('-r', '--reverse', action='store_true', help='reverse (family->superfamily)')
	parser.add_argument('infile', help='plaintext file to convert')
	parser.add_argument('outfile', help='where to write the resulting json')

	args = parser.parse_args()

	main(args.infile, args.outfile, reverse=args.reverse)
