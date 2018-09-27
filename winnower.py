#!/usr/bin/env python

from __future__ import print_function, division, generators

import argparse, os, sys, json
import extendomatic

class Alignment(object):
	def __init__(self, name, jstr, tmdatadir=None, dthreshold=4):
		self.name = name
		self.query, self.qchain, self.qhel, vs, self.subject, self.schain, self.shel = name.split('_')
		self.qhels = get_hel_list(self.qhel)
		self.shels = get_hel_list(self.shel)
		obj = json.loads(jstr)
		self.jstr = jstr.strip()
		data = extendomatic.unpack_obj(obj, dthreshold=dthreshold)
		self.rmsd = data[0]
		self.length = data[1]
		self.mincov = data[2]
		self.maxcov = data[3]

		self.qpresent = obj['qpresent']
		self.spresent = obj['spresent']

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

def load_opm_tmss(tmdatadir='tmdata'):
	tmdict = {}
	with open('{}/opm/ASSIGNMENTS.TSV'.format(tmdatadir)) as f:
		for l in f: 
			pdbid, tmstr = l.strip().split('\t')
			pdbid = pdbid[:4].upper() + pdbid[4:]
			tmdict[pdbid] = []
			tmdict[pdbid[:4]] = []
			for span in tmstr.split(','):
				tmdict[pdbid].append([int(x) for x in span.split('-')])
				tmdict[pdbid[:4]].append([int(x) for x in span.split('-')])
	return tmdict

def overlap(s1, s2):
	r1 = set(range(*s1))
	r2 = set(range(*s2))
	return len(r1.intersection(r2))
def count_covered_tmss(spans, interval, threshold=10):
	tmss = 0
	for span in spans:
		n = overlap(span, interval)
		if n > threshold: tmss += 1
	return tmss

def main(infile, outfile='/dev/stdout', stretch=0, append=False, tmdatadir='tmdata', dthreshold=4, mincov=0.):
	print(infile, outfile)
	if not tmdatadir: tmdatadir = None
	currstrucs = (None, None)

	tmdict = load_opm_tmss(tmdatadir)

	donehels = ([], [])
	#best = {}
	best = []
	if append: g = open(outfile, 'a')
	else: g = open(outfile, 'w')
	with open(infile) as f:
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue

			try: name, jstr = l.split('\t')
			except ValueError:
				print(l)
				continue
			try: current = Alignment(name, jstr, tmdatadir=tmdatadir, dthreshold=dthreshold)
			except ValueError: continue
			if current.rmsd == -1: continue
			if current.mincov < mincov: continue
			#print(name, current.rmsd, current.mincov)
			
			query, qchain, qhel, vs, subject, schain, shel = name.split('_')
			qhels = get_hel_list(qhel)
			shels = get_hel_list(qhel)
			qpdbid = '{}_{}'.format(query, qchain)
			spdbid = '{}_{}'.format(subject, schain)

			#eliminate alignments that couldn't ever exceed 2 TMSs
			skip = False
			try: qtmscount = count_covered_tmss(tmdict[qpdbid], current.qpresent[0])
			except KeyError: 
				try: qtmscount = count_covered_tmss(tmdict[qpdbid[:4]], current.qpresent[0])
				except KeyError: qtmscount = 6.5 #arbitrary value to allow PDBTM-only structures for now
			finally:
				if qtmscount <= 2: skip = True
			try: stmscount = count_covered_tmss(tmdict[spdbid], current.spresent[0])
			except KeyError: 
				try: stmscount = count_covered_tmss(tmdict[spdbid[:4]], current.spresent[0])
				except KeyError: stmscount = 6.5 #arbitrary value to allow PDBTM-only structures for now
			finally:
				if stmscount <= 2: skip = True
			if skip: continue

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
				currstrucs = (qpdbid, spdbid)

				best.append(current)
			else:
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
	parser.add_argument('-d', default=4, type=float, help='distance threshold (default:4)')
	parser.add_argument('--min-cov', default=0., type=float, help='minimum coverage')
	parser.add_argument('-s', default=0, type=int, help='how many extra alignments to keep for each pair of structures (default:0)')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout')
	parser.add_argument('--tmdatadir', nargs='?', default='tmdata')

	args = parser.parse_args()
	if args.infile == args.outfile: 
		parser.print_usage()
		exit(1)

	if not os.path.isdir(args.tmdatadir): 
		raise IOError('Could not find tmdata directory')

	if os.path.isdir(args.infile):
		open(args.outfile, 'w')
		for famvfam in os.listdir(args.infile):
			for spfn in os.listdir('{}/{}/superpositions'.format(args.infile, famvfam)):
				if spfn.startswith('.'): continue
				elif not spfn.lower().endswith('tsv'): continue
				main('{}/{}/superpositions/{}'.format(args.infile, famvfam, spfn), args.outfile, stretch=args.s, append=True, tmdatadir=args.tmdatadir, dthreshold=args.d, mincov=args.min_cov)
	else: main(args.infile, args.outfile, stretch=args.s, append=False, tmdatadir=args.tmdatadir, dthreshold=args.d, mincov=args.min_cov)
