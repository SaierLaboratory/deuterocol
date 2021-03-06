#!/usr/bin/env python

from __future__ import print_function, division, generators

import argparse, os, sys, json
import extendomatic

class Alignment(object):
	def __init__(self, name, jstr, tmdatadir=None, dthreshold=4, scoring_function='mincov'):
		self.name = name
		self.query, self.qchain, self.qhel, vs, self.subject, self.schain, self.shel = name.split('_')
		self.qhels = get_hel_list(self.qhel)
		self.shels = get_hel_list(self.shel)
		obj = json.loads(jstr)
		self.jstr = jstr.strip()
		data = extendomatic.unpack_obj(obj, dthreshold=dthreshold)
		self.rmsd = data['rmsd']
		self.length = data['length']
		self.mincov = data['mincov']
		self.maxcov = data['maxcov']
		self.quality = obj['quality']

		self.qpresent = obj['qpresent']
		self.spresent = obj['spresent']

		self.scoring_function = scoring_function

	#def __lt__(self, other):
	#	''' maximizes coverage '''
	#	if self.mincov > other.mincov: return True
	#	elif self.mincov == other.mincov:
	#		if self.rmsd < other.rmsd: return True
	#		else: return False
	#	else: return False

	def _score_mincov(self): 
		''' maximizes coverage '''
		return self.mincov

	def _score_cartesian(self, covweight=1): 
		''' minimizes distance to (RMSD=0.0, mincov=1.0) '''
		return (self.rmsd)**2 + covweight*(1.0 - self.mincov)**2

	def _score_quality(self):
		''' maximizes quality '''
		return -self.quality

	def get_score(self):
		''' scores an alignment '''
		#return self._score_cartesian(covweight=8)
		#return self._score_quality()
		if self.scoring_function.startswith('mincov'): return self._score_mincov()
		elif self.scoring_function.startswith('cartes'): return self._score_cartesian(covweight=8)
		elif self.scoring_function.startswith('qual'): return self._score_quality()
		else: raise ValueError('Invalid scoring function {}'.format(self.scoring_function))

	def __lt__(self, other):
		if self.get_score() < other.get_score(): return True
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

def main(infile, outfile='/dev/stdout', stretch=0, append=False, tmdatadir='tmdata', dthreshold=4, mincov=0., scoring_function='mincov'):
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
			try: current = Alignment(name, jstr, tmdatadir=tmdatadir, dthreshold=dthreshold, scoring_function=scoring_function)
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
	parser.add_argument('-d', default=6, type=float, help='distance threshold for coverage calculations (float, default:6)')
	parser.add_argument('--min-cov', default=0., type=float, help='minimum coverage (float [0..1], default:0.)')
	parser.add_argument('-s', default=0, type=int, help='how many extra alignments to keep for each pair of structures (int, default:0)')
	parser.add_argument('outfile', nargs='?', default='/dev/stdout', help='where to write the surviving alignments')
	parser.add_argument('--tmdatadir', nargs='?', default='tmdata', help='where the TMdata assignments are stored')
	parser.add_argument('--scoring-function', default=None, help='scoring function to use (options: qual, mincov, cartes)')

	args = parser.parse_args()
	if args.infile == args.outfile: 
		parser.print_usage()
		exit(1)

	if args.scoring_function is None:
		parser.print_usage()
		exit(1)
	elif args.scoring_function not in ['mincov', 'cartes', 'qual']: 
		parser.print_usage()
		exit(1)
		

	if not os.path.isdir(args.tmdatadir): 
		raise IOError('Could not find tmdata directory')

	#TODO: expose this to argparser
	subpath = 'tmalignments'

	if os.path.isdir(args.infile):
		open(args.outfile, 'w')
		for famvfam in os.listdir(args.infile):
			if '_vs_' not in famvfam: continue
			if not os.path.isdir('{}/{}/{}'.format(args.infile, famvfam, subpath)): 
				print('[WARNING]: Could not find TMalign results in directory {}/{}/{}'.format(args.infile, famvfam, subpath), file=sys.stderr)
				continue
			for spfn in os.listdir('{}/{}/{}'.format(args.infile, famvfam, subpath)):
				if spfn.startswith('.'): continue
				elif not spfn.lower().endswith('tsv'): continue
				main('{}/{}/{}/{}'.format(args.infile, famvfam, subpath, spfn), args.outfile, stretch=args.s, append=True, tmdatadir=args.tmdatadir, dthreshold=args.d, mincov=args.min_cov, scoring_function=args.scoring_function)
	else: main(args.infile, args.outfile, stretch=args.s, append=False, tmdatadir=args.tmdatadir, dthreshold=args.d, mincov=args.min_cov, scoring_function=args.scoring_function)
