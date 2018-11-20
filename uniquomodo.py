#!/usr/bin/env python

from __future__ import print_function, division
import sys
import json
import argparse
import shutil
import tempfile

import pymol

def test_transform():

	matrix = [
     0.9164623314, -0.2389936591,  0.3209031413, 495.2727054608,  
     0.1849898163, -0.4580664883, -0.8694560714, 132.6568206768,  
     0.3547894629,  0.8561875514, -0.3755892888, 186.3692863698,  
     0.000000000000,  0.0000000000,  0.0000000000,  1.0000000000, 
	]

	#matrix = [
	#	10., 0., 0., 0.,
	#	0., 1., 0., 0.,
	#	0., 0., 1., 0.,
	#	0., 0., 0., 1.,
	#]

	pymol.finish_launching()
	pymol.cmd.load('2v8n.pdb')
	pymol.cmd.load('4pyp.pdb')
	help(pymol.cmd.transform_selection)
	pymol.cmd.transform_selection('2v8n', matrix)
	pymol.cmd.center('2v8n')
	pass


def extract_selectors(name, obj):
	query, qchain, qtms, vs, subject, schain, stms = name.split('_')
	seldict = {}

	seldict['qchain'] = '{} and c. {}'.format(query, qchain)
	seldict['schain'] = '{} and c. {}'.format(subject, schain)
	seldict['qmasked'] = '{} and c. {} and i. {}-{}'.format(
		obj['query'], 
		qchain,
		obj['qpresent'][0][0],
		obj['qpresent'][0][1]
	)
	seldict['smasked'] = '{} and c. {} and i. {}-{}'.format(
		obj['subject'], 
		schain,
		obj['spresent'][0][0],
		obj['spresent'][0][1]
	)
	alnq = []
	for span in obj['qaligned']:
		if span[0] is None: alnq += [None] * span[1]
		else: alnq += list(range(span[0], span[0]+span[1]+1))
	alns = []
	for span in obj['saligned']:
		if span[0] is None: alns += [None] * span[1]
		else: alns += list(range(span[0], span[0]+span[1]+1))

	seldict['qpresent'] = '{} and c. {} and i. '.format(obj['query'], qchain)
	seldict['spresent'] = '{} and c. {} and i. '.format(obj['subject'], schain)
	seldict['qaligned'] = '{} and c. {} and i. '.format(obj['query'], qchain)
	seldict['saligned'] = '{} and c. {} and i. '.format(obj['subject'], schain)

	qp = obj['qpresent'][0][1] - obj['qpresent'][0][0]
	sp = obj['spresent'][0][1] - obj['spresent'][0][0]
	qa = 0
	sa = 0

	for q, d, s in zip(alnq, obj['distances'], alns):
		if q is not None: seldict['qpresent'] += '{}+'.format(q)
		if s is not None: seldict['spresent'] += '{}+'.format(s)
		if (q is not None) and (s is not None) and (d is not None):
			seldict['qaligned'] += '{}+'.format(q)
			seldict['saligned'] += '{}+'.format(s)
			qa += 1
			sa += 1

	print('='*80)
	print('{}:'.format(name))
	print('\tRMSD: {:0.2f}'.format(obj['rmsd']))
	if qa/qp < sa/sp:
		print('\tMinimum coverage: {:0.2%} of query'.format(qa/qp))
	else:
		print('\tMinimum coverage: {:0.2%} of subject'.format(sa/sp))
	seldict['qpresent'] = seldict['qpresent'][:-1]
	seldict['spresent'] = seldict['spresent'][:-1]
	seldict['qaligned'] = seldict['qaligned'][:-1]
	seldict['saligned'] = seldict['saligned'][:-1]
	return seldict

def main(name, obj, wd=None):
	sys.argv = sys.argv[:1]
	pymol.finish_launching()

	if wd is None:
		print('Error: No output directory specified')
		exit(1)

	pymol.cmd.set('fetch_path', wd)
	pymol.cmd.fetch(obj['query'])
	pymol.cmd.fetch(obj['subject'])

	selectors = extract_selectors(name, obj)
	pymol.cmd.hide('lines')
	for k in sorted(selectors):
		if k != 'matrix':
			pymol.cmd.select(k, selectors[k])
			if k.endswith('masked'):
				pymol.cmd.show('cartoon', k)

	truematrix = []
	for row in obj['matrix']: truematrix += row
	truematrix += [0., 0., 0., 1.]
	pymol.cmd.transform_selection(obj['query'], truematrix)
	pymol.cmd.rotate('x', 90)
	pymol.cmd.center('*masked')
	pymol.cmd.zoom('*masked')

	pymol.cmd.color('palegreen', obj['query'])
	pymol.cmd.color('palecyan', obj['subject'])
	pymol.cmd.color('tv_green', 'qaligned')
	pymol.cmd.color('cyan', 'saligned')

	pymol.cmd.deselect()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', nargs='+', help='TSVs to read data from')

	parser.add_argument('name', help='name of alignment to load. Be warned that this selects only the first occurrence in all the files in the order given')

	args = parser.parse_args()

	found = False
	jstr = ''
	for fn in args.infile:
		with open(fn) as f:
			for l in f:
				if l.startswith(args.name):
					print(args.name, 'found in', fn)
					jstr = l[l.find('\t')+1:]
				if found: break
		if found: break

	if not jstr: 
		print('Error: Could not find {} in {}'.format(args.name, args.infile))

	obj = json.loads(jstr)
	wd = tempfile.mkdtemp()
	try:
		main(args.name, obj, wd=wd)
	finally:
		shutil.rmtree(wd)
