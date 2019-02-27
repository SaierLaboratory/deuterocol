#!/usr/bin/env python
from __future__ import print_function, division, unicode_literals

import re
import argparse
import sys
import os

try: from Bio.PDB import protein_letters_3to1 as CODE
except AttributeError: from Bio.PDB import to_one_letter_codes as CODE

class NumberedSequence(object):
	def __init__(self, iterable):
		self.iterable = iterable
		self.validate()

	def validate(self):
		for i, pair in enumerate(self.iterable):
			if not ((len(pair) == 2) and (type(pair[1]) is int)): 
				raise ValueError('Element {} is invalid: {}'.format(i, pair))


	def get_by_resi(self, start, end=None, step=1): 
		if end is None: index = range(start, start+1)
		else: index = range(start, end+1, step)

		seq = ''
		for pair in self:
			if pair[1] in index: seq += pair[0]
		return seq
		#else: raise IndexError('list index out of range')

	def get_seq(self): return ''.join([pair[0] for pair in self])

	def get_contigs(self):
		rawcontigs = []
		lasti = None
		for pair in self:
			if lasti is None: rawcontigs.append([pair])
			elif pair[1] == (lasti + 1): rawcontigs[-1].append(pair)
			else: rawcontigs.append([pair])
			lasti = pair[1]
		return [NumberedSequence(contig) for contig in rawcontigs]

	def __getitem__(self, index):
		if type(index) is int:
			for pair in self:
				if pair[1] == index: return NumberedSequence([self.iterable[index]])
		else:
			start = self.iterable[0][1] if index.start is None else index.start
			stop = self.iterable[-1][1] if index.stop is None else index.stop
			step = 1 if index.step is None else index.step
			requested = range(start, stop+step, step)
			substring = []
			for pair in self:
				if pair[1] in requested: substring.append(pair)
			return NumberedSequence(substring)
	def find(self, pattern, start=None):
		def _sift(self, pattern, bank=None, start=None):
			if bank is None:
				bank = []
				for pair in self:
					if (start) and (pair[1] < start): continue
					if pair[0] == pattern[0]: 
						bank.append([pair])
				return bank
			else: 
				bank = bank[:]
				deleteme = []
				for segnum, segment in enumerate(bank):
					resi = segment[-1][1]+1
					resn = self.get_by_resi(resi)
					if resn != pattern: deleteme.append(segnum)
					else: segment.append((resn, resi))

				for segnum in deleteme[::-1]: bank.pop(segnum)
				return bank
		bank = None
		for length in range(len(pattern)):
			bank = _sift(self, pattern[length], bank, start=start)
		if bank: return NumberedSequence(bank[0])
		else: return None

	def gapless_align_string(self, string):
		if len(string) > len(self): raise NotImplementedError
		#elif len(string) == len(self):
		#	idents = 0
		#	for c1, c2 in zip(self, string): 
		#		if c1 == c2: idents += 1
		#	if idents >= 0.9 * len(self): return self[:]
		else:
			idents = {}
			maxident = 0
			#print(self, string)
			for i in range(len(self) - len(string) + 1):
				idents[i] = 0
				for c1, c2 in zip(self.iterable[i:i+len(string)], string):
					if c1[0] == c2: idents[i] += 1
				maxident = max(maxident, idents[i])
				if maxident == len(string): 
					return NumberedSequence(self.iterable[i:i+len(string)])
			#print('='*80)
			#print(self, string)
			#for i in idents: print(i, idents[i])
			for i in idents:
				if idents[i] == maxident: return NumberedSequence(self.iterable[i:i+len(string)])

	def get_range(self):
		start, stop = None, None
		for pair in self:
			if start is None: start = pair[1]
			stop = pair[1]
		return start, stop
		

	def __iter__(self): return iter(self.iterable)
	def __repr__(self):
		low, high = None, None
		seq = ''
		for pair in self:
			if low is None: low = pair[1]
			seq += pair[0]
			high = pair[1]
		return '<NumberedSequence range=({}, {}) seq="{}">'.format(low, high, seq)
	def __len__(self): return len(self.iterable)

def extract_atom_sequence(f, chain):
	sequence = []
	for l in f:
		if l.startswith('ATOM') or l.startswith('HETATM'):
			if 'CA' in l: 
				if l[21] != chain: continue
				try: resn = CODE[l[17:20]]
				except KeyError: continue
				try: resi = int(l[22:26])
				except ValueError: print(l)
				sequence.append((resn, resi))
	return NumberedSequence(sequence)

def get_aligned_contig_sequences(alignmentlines, minlength=4):
	lastmidline = None
	qcontigs = []
	tcontigs = []
	permitted = ':'
	for qresn, midline, tresn in zip(*alignmentlines):
		if lastmidline is None and midline in ' ': pass
		elif lastmidline is None and midline in permitted:
			qcontigs.append(qresn)
			tcontigs.append(tresn)
		elif lastmidline in permitted and midline in permitted:
			qcontigs[-1] += qresn
			tcontigs[-1] += tresn
		elif lastmidline in ' ' and midline in permitted:
			qcontigs.append(qresn)
			tcontigs.append(tresn)
		lastmidline = midline
	filtered_qcontigs = []
	filtered_tcontigs = []
	for qcontig in qcontigs:
		if len(qcontig) >= minlength: filtered_qcontigs.append(qcontig)
	for tcontig in tcontigs:
		if len(tcontig) >= minlength: filtered_tcontigs.append(tcontig)
	return filtered_qcontigs, filtered_tcontigs

def main(rawargs):
	parser = argparse.ArgumentParser()
	parser.add_argument('infile')
	args = parser.parse_args(rawargs)

	fns = []
	alignmentlines = []
	recording = False
	with open(args.infile) as f:
		for l in f:
			if not l.strip(): continue
			elif l.startswith('Name of'):
				fns.append(re.split('\s*:\s*', l.strip())[1])
			elif l.startswith('(":" denotes'): recording = True
			elif recording: alignmentlines.append(l.replace('\n', ''))

	pdbsequences = []
	chains = []
	pdbaccs = []
	for fn in fns: 
		bn = os.path.basename(fn)
		pdbaccs.append(bn[:4])
		chain = re.findall('(.)_h', bn)[0]
		chains.append(chain)
		pdbsequences.append(extract_atom_sequence(open(fn), chain=chain))

	qcontigs, tcontigs = get_aligned_contig_sequences(alignmentlines, minlength=4)
	rangestr = ''

	selectors = []
	for qcontig in qcontigs:
		#contobj = pdbsequences[0].find(qcontig.replace('-', ''))
		contobj = pdbsequences[0].gapless_align_string(qcontig.replace('-', ''))
		rangestr += '{}-{}+'.format(*contobj.get_range())
	rangestr = rangestr[:-1]
	selectors.append('{} and c. {} and i. {}'.format(pdbaccs[0], chains[0], rangestr))
		
	rangestr = ''
	for tcontig in tcontigs:
		#contobj = pdbsequences[1].find(tcontig.replace('-', ''))
		contobj = pdbsequences[1].gapless_align_string(tcontig.replace('-', ''))
		rangestr += '{}-{}+'.format(*contobj.get_range())
	rangestr = rangestr[:-1]
	selectors.append('{} and c. {} and i. {}'.format(pdbaccs[1], chains[1], rangestr))
	return selectors
	#for i in range(0, len(alignmentlines[0]), 80):
	#	print('q', alignmentlines[0][i:i+80])
	#	print('m', alignmentlines[1][i:i+80])
	#	print('t', alignmentlines[2][i:i+80])

if __name__ == '__main__':
	main(sys.argv[1:])
