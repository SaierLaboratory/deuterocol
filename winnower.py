#!/usr/bin/env python

from __future__ import print_function, division, generators

import argparse, os, sys, json
import extendomatic

class Alignment(object):
	def __init__(self, name, jstr):
		self.name = name
		self.query, self.qchain, self.qhel, vs, self.subject, self.schain, self.shel = name.split('_')
		self.qhels = get_hel_list(self.qhel)
		self.shels = get_hel_list(self.shel)
		obj = json.loads(jstr)
		self.jstr = jstr.strip()
		self.rmsd, self.length, self.mincov, self.maxcov = extendomatic.unpack_obj(obj)
	def __lt__(self, other):
		if self.mincov > other.mincov: return True
		elif self.mincov == other.mincov:
			if self.rmsd < other.rmsd: return True
			else: return False
		else: return False

def get_hel_list(s):
	if s.startswith('h'): s = s[1:]
	start, end = s.split('-')
	return list(range(int(start), int(end)+1))

def intersects(l1, l2):
	for x1 in l1:
		if x1 in l2: return True
	return False


def main(infile, outfile='/dev/stdout', stretch=0, append=False):
	currstrucs = (None, None)
	donehels = ([], [])
	#best = {}
	best = []
	if append: g = open(outfile, 'a')
	else: g = open(outfile, 'w')
	with open(infile) as f:
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue

			name, jstr = l.split('\t')
			try: current = Alignment(name, jstr)
			except ValueError: continue
			query, qchain, qhel, vs, subject, schain, shel = name.split('_')
			qhels = get_hel_list(qhel)
			shels = get_hel_list(qhel)
			qpdbid = '{}_{}'.format(query, qchain)
			spdbid = '{}_{}'.format(subject, schain)
			
		#	if currstrucs == (qpdbid, spdbid):
		#		if best[currstrucs][0].mincov < current.mincov:
		#			best[currstrucs].insert(0, current)
		#		elif best[currstrucs][0].mincov == current.mincov:
		#			if best[currstrucs][0].rmsd > current.rmsd:
		#				best[currstrucs].insert(0, current)
		#			else:
		#				best[currstrucs].insert(1, current)
		#		else: best[currstrucs].append(current)
		#		if len(best[currstrucs]) > (stretch + 1): best[currstrucs].pop(-1)
		#	else:
		#		if currstrucs[0]: 
		#			for aln in best[currstrucs]:
		#				print('{}\t{}'.format(aln.name, aln.jstr), file=g)
		#		currstrucs = (qpdbid, spdbid)
		#		best[currstrucs] = [current]
			if currstrucs == (qpdbid, spdbid):
				print('if currstrucs == (qpdbid, spdbid)')
				#if best[0].mincov < current.mincov:
				#	best.insert(0, current)
				#elif best[0].mincov == current.mincov:
				#	if best[0].rmsd > current.rmsd:
				#		best.insert(0, current)
				#	else:
				#		best.insert(1, current)
				#else: best.append(current)
				#if len(best) > (stretch + 1): best.pop(-1)
				best.append(current)
			elif currstrucs[0] is None:
				print('elif currstrucs[0] is None')
				currstrucs = (qpdbid, spdbid)

				best.append(current)
			else:
				print('else')
				donehels = [[], []]
				for aln in sorted(best)[:stretch+1]:
					if intersects(aln.qhels, donehels[0]) or intersects(aln.shels, donehels[1]): continue
					print('{}\t{}'.format(aln.name, aln.jstr), file=g)
					donehels[0] += aln.qhels
					donehels[1] += aln.shels
				#print(currstrucs)
				best = [current]
				currstrucs = (qpdbid, spdbid)
				donehels = (qhels, shels)
				


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', help='a single superposition TSV or a deuterocol2 root directory')
	parser.add_argument('-s', default=0, type=int, help='how many alignments to keep for each pair of structures')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout')

	args = parser.parse_args()
	if args.infile == args.outfile: 
		parser.print_usage()
		exit(1)

	if os.path.isdir(args.infile):
		open(args.outfile, 'w')
		for famvfam in os.listdir(args.infile):
			for spfn in os.listdir('{}/{}/superpositions'.format(args.infile, famvfam)):
				if spfn.startswith('.'): continue
				elif not spfn.lower().endswith('tsv'): continue
				main('{}/{}/superpositions/{}'.format(args.infile, famvfam, spfn), args.outfile, stretch=args.s, append=True)
	else: main(args.infile, args.outfile, stretch=args.s, append=False)
