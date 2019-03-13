#!/usr/bin/env python

from __future__ import print_function, division
import sys
import os
import json
import argparse
import shutil
import tempfile
import subprocess

import pymol

import tmalignparser

def error(*stuff):
	print('[ERROR]:', *stuff, file=sys.stderr)
	sys.exit(1)

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


def to_pymol_ranges(spans):
	out = ''
	for span in spans:
		if span[0] == span[1]:  out += '{}+'.format(span[0])
		else: out += '{}-{}+'.format(span[0], span[-1])
	return out
		
def compress_list(l):
	lasti = None
	spans = []
	for i in l:
		if lasti is None: spans.append([i, i])
		elif i == (lasti + 1): spans[-1][1] = i
		elif i == lasti: pass
		else: spans.append([i, i])
		lasti = i

	return spans

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

	#this will get ignored anyway
	seldict['qaligned'] = '{} and c. {} and i. '.format(obj['query'], qchain)
	seldict['saligned'] = '{} and c. {} and i. '.format(obj['subject'], schain)

	qp = obj['qpresent'][0][1] - obj['qpresent'][0][0]
	sp = obj['spresent'][0][1] - obj['spresent'][0][0]
	qa = 0
	sa = 0

	qpresent = []
	spresent = []

	for q, d, s in zip(alnq, obj['distances'], alns):
		if q is not None: qpresent.append(q)
		if s is not None: spresent.append(s)
		if (q is not None) and (s is not None) and (d is not None):
			qpresent.append(q)
			spresent.append(s)
			qa += 1
			sa += 1

	seldict['qpresent'] += to_pymol_ranges(compress_list(qpresent))
	seldict['spresent'] += to_pymol_ranges(compress_list(spresent))

	print('='*80)
	print('{}:'.format(name))
	print('\tRMSD: {:0.2f}'.format(obj['rmsd']))
	print('\tTM-score: {:0.2f}'.format(obj['quality']))
	if qa/qp < sa/sp:
		print('\tMinimum coverage: {:0.2%} of query'.format(qa/qp))
		print('\tMaximum coverage: {:0.2%} of subject'.format(sa/sp))
	else:
		print('\tMinimum coverage: {:0.2%} of subject'.format(sa/sp))
		print('\tMaximum coverage: {:0.2%} of query'.format(qa/qp))
	seldict['qpresent'] = seldict['qpresent'][:-1]
	seldict['spresent'] = seldict['spresent'][:-1]
	seldict['qaligned'] = seldict['qaligned'][:-1]
	seldict['saligned'] = seldict['saligned'][:-1]
	return seldict

def run_tmap(obj, d2dir):
	qfn = '{}/cut_pdbs/{}'.format(d2dir, os.path.basename(obj['queryfn']))
	sfn = '{}/cut_pdbs/{}'.format(d2dir, os.path.basename(obj['subjectfn']))

	tf = tempfile.NamedTemporaryFile()
	out = subprocess.check_output(['TMalign', qfn, sfn])
	tf.write(out)
	tf.flush()

	selectors = tmalignparser.main([tf.name, qfn, sfn])
	return selectors

def main(name, obj, d2dir, wd=None):
	sys.argv = sys.argv[:1]
	pymol.finish_launching()

	if wd is None:
		print('Error: No output directory specified')
		exit(1)

	pymol.cmd.set('fetch_path', wd)
	#pymol.cmd.fetch(obj['query'])
	pymol.cmd.load('{}/pdbs/{}.pdb'.format(d2dir, obj['query']))
	#pymol.cmd.fetch(obj['subject'])
	pymol.cmd.load('{}/pdbs/{}.pdb'.format(d2dir, obj['subject']))

	selectors = extract_selectors(name, obj)

	qaligned, saligned = run_tmap(obj, d2dir)
	selectors['qaligned'] = qaligned
	selectors['saligned'] = saligned

	pymol.cmd.hide('lines')
	print()
	print('Selectors:')
	print('----------')
	for k in sorted(selectors):
		if k != 'matrix':
			pymol.cmd.select(k, selectors[k])
			
			if k.endswith('masked'): 
				pymol.cmd.show('cartoon', k)
			print(k + ':', selectors[k])


	truematrix = []
	for row in obj['matrix']: truematrix += row
	truematrix += [0., 0., 0., 1.]
	pymol.cmd.transform_selection(obj['query'], truematrix)
	pymol.cmd.rotate('x', 90)
	pymol.cmd.center('*present')
	pymol.cmd.zoom('*present')
	pymol.cmd.show('cartoon', '*present')

	pymol.cmd.color('palegreen', obj['query'])
	pymol.cmd.color('palecyan', obj['subject'])
	pymol.cmd.color('tv_green', 'qaligned')
	pymol.cmd.color('cyan', 'saligned')

	pymol.cmd.deselect()

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('d1dir', help='Deuterocol1 directory')
	parser.add_argument('d2dir', help='Deuterocol2 directory')

	parser.add_argument('name', help='name of alignment to load. Be warned that this selects only the first occurrence in all the files in the order given')

	args = parser.parse_args()

	query, qchain, qhel, vs, subject, schain, shel = args.name.split('_')
	possfams = set()

	with open(args.d1dir + '/pdblist.json') as f: pdbmap = json.loads(f.read())
	for fam in pdbmap:
		if (query + '_' + qchain) in pdbmap[fam]: possfams.add(fam)
		if (subject + '_' + schain) in pdbmap[fam]: possfams.add(fam)

	possfams = tuple(possfams)
	subdirlist = []
	for subdir in os.listdir(args.d2dir):
		count = 0
		for fam in possfams: 
			if fam in subdir: count += 1
		if count >= 2: subdirlist.append(subdir)
	if not subdirlist: error('Could not find directories for comparisons between families {}'.format(possfams))

	found = False
	jstr = ''
	n = 0

	for subdir in subdirlist:
		fn = args.d2dir + '/' + subdir + '/tmalignments/sp_all.tsv'
		if not os.path.isfile(fn): error('Could not find TMalign records for', subdir)

	#for fn in args.infile:
		with open(fn) as f:
			for l in f:
				n += 1
				if not n % 10000: print('Checked', n, 'lines')
				if l.startswith(args.name[:2]):
					if l.startswith(args.name):
						print('Found', args.name, 'in', fn)
						jstr = l[l.find('\t')+1:]
						found = True
				if found: break
		if found: break

	if not jstr: 
		print('Error: Could not find {} in {}'.format(args.name, args.infile))

	obj = json.loads(jstr)
	wd = tempfile.mkdtemp()
	try: main(args.name, obj, wd=wd, d2dir=args.d2dir)
	finally: shutil.rmtree(wd)
