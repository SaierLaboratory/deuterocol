#!/usr/bin/env python2
from __future__ import division

import sys, json
import numpy as np

def get_coverages(obj):
	qaligned = 0
	for span in obj['qaligned']:
		if span[0] is None: continue
		else: qaligned += span[1] - span[0] + 1
	saligned = 0
	for span in obj['saligned']:
		if span[0] is None: continue
		else: saligned += span[1] - span[0] + 1
	qpresent = obj['qpresent'][0][1] - obj['qpresent'][0][0] + 1
	spresent = obj['spresent'][0][1] - obj['spresent'][0][0] + 1
	return qaligned/qpresent, saligned/spresent

rmsds = []
lengths = []
qcovs = []
scovs = []
mincovs = []
with open(sys.argv[1]) as f:
	for l in f:
		sl = l.split('\t')
		try: obj = json.loads(sl[1])
		except ValueError: continue
		#everything: R^2 is 0.09
		if obj['rmsd'] == -1: continue
		if obj['rmsd'] <= 4.5 and obj['length'] >= 60:
			rmsds.append(obj['rmsd'])
			lengths.append(obj['length'])
			#print('{}\t{}'.format(obj['rmsd'], obj['length']))

			qcov, scov = get_coverages(obj)
			qcovs.append(qcov)
			scovs.append(scov)
			mincovs.append(min(qcov, scov))

def diffscore(x, y):
	X = np.array(x)
	Y = np.array(y)
	ZX = (X - np.mean(X))/np.std(X)
	ZY = (Y - np.mean(Y))/np.std(Y)
	Z = ZY**2*np.sign(ZY) - ZX**2*np.sign(ZX)

	return Z

def newscoringfunction(x, y):
	X = np.array(x)
	Y = np.array(y)

	Z = 1
	return Z

scores = diffscore(rmsds, mincovs)

for rmsd, mincov, score in zip(rmsds, mincovs, scores):
	print('{}\t{:0.4f}\t{:0.2f}'.format(rmsd, mincov, score))
