#!/usr/bin/env python

from __future__ import print_function, division
import json, os
import tempfile, urllib
import subprocess
import pymol
import argparse
import sys
def extract_selectors(name, obj):
	query, qchain, qtms, vs, subject, schain, stms = name.split('_')

	qmasked = '{} and c. {} and i. {}-{}'.format(
		obj['query'],
		qchain,
		obj['qpresent'][0][0],
		obj['qpresent'][0][1],
	)
	smasked = '{} and c. {} and i. {}-{}'.format(
		obj['subject'],
		schain,
		obj['spresent'][0][0],
		obj['spresent'][0][1],
	)

	alnq = []
	for span in obj['qaligned']:
		if span[0] is None: alnq += [None] * span[1]
		else: alnq += list(range(span[0], span[1]+1))
	alns = []
	for span in obj['saligned']:
		if span[0] is None: alns += [None] * span[1]
		else: alns += list(range(span[0], span[1]+1))

	qpresent = '{} and c. {} and i. '.format(obj['query'], qchain)
	spresent = '{} and c. {} and i. '.format(obj['subject'], schain)
	qaligned = '{} and c. {} and i. '.format(obj['query'], qchain)
	saligned = '{} and c. {} and i. '.format(obj['subject'], schain)

	qp = obj['qpresent'][0][1] - obj['qpresent'][0][0]
	sp = obj['spresent'][0][1] - obj['spresent'][0][0]
	qa = 0
	sa = 0
	print(obj)
	print(obj['distances'])
	print(alnq, alns)
	for q, d, s in zip(alnq, obj['distances'], alns):
		if q is not None: qpresent += '{}+'.format(q)
		if s is not None: spresent += '{}+'.format(s)
		if (q is not None) and (s is not None) and (d is not None):
			qaligned += '{}+'.format(q)
			saligned += '{}+'.format(s)
			qa += 1
			sa += 1

	print('='*80)
	print('{}:'.format(name))
	print('\tRMSD: {:0.2f}'.format(obj['rmsd']))
	if qa/qp < sa/sp:
		print('\tMinimum coverage: {:0.2%} of query'.format(qa/qp))
	else:
		print('\tMinimum coverage: {:0.2%} of subject'.format(sa/sp))
	qpresent = qpresent[:-1]
	spresent = spresent[:-1]
	qaligned = qaligned[:-1]
	saligned = saligned[:-1]

	return qmasked, smasked, qpresent, spresent, qaligned, saligned

def main(name, infile):
	with open(infile) as f:
		found = False
		for l in f:
			if l.startswith(name): 
				found = l.strip()
		if not found:
			print('Error: Could not find {} in {}'.format(name, infile))
			exit(1)
	name, jstr = found.split('\t')
	obj = json.loads(jstr)

	prefix = tempfile.mkdtemp()

	query, qchain, qtms, vs, subject, schain, stms = name.split('_')

	f = urllib.urlopen('https://files.rcsb.org/view/{}.pdb'.format(obj['query']))
	with open('{}/{}_preimage.pdb'.format(prefix, obj['query']), 'w') as g: g.write(f.read())
	f = urllib.urlopen('https://files.rcsb.org/view/{}.pdb'.format(obj['subject']))
	with open('{}/{}.pdb'.format(prefix, obj['subject']), 'w') as g: g.write(f.read())
	p = subprocess.Popen(['superpose', 
		'{}/{}_preimage.pdb'.format(prefix, obj['query']), 
		'-s', '{}/{}-{}'.format(qchain, obj['qpresent'][0][0], obj['qpresent'][0][1]), 
		'{}/{}.pdb'.format(prefix, obj['subject']), 
		'-s', '{}/{}-{}'.format(schain, obj['spresent'][0][0], obj['spresent'][0][1]),
		'-o', '{}/{}.pdb'.format(prefix, obj['query'])
	], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()

	for l in out:
		if 'r.m.s.d' in l: 
			rmsd = l.strip().split()[1]
			if rmsd != obj['rmsd']:
				print('[WARN]: Rerun RMSD ({}) does not match recorded RMSD ({})')
			break

	sys.argv = sys.argv[:1]

	pymol.finish_launching()
	pymol.cmd.load('{}/{}.pdb'.format(prefix, obj['query']))
	pymol.cmd.load('{}/{}.pdb'.format(prefix, obj['subject']))

	qmasked, smasked, qpresent, spresent, qaligned, saligned = extract_selectors(name, obj)

	def spview():
		pymol.cmd.select('qchain', '{} and c. {}'.format(query, qchain))
		pymol.cmd.select('schain', '{} and c. {}'.format(subject, schain))
		pymol.cmd.select('qmasked', qmasked)
		pymol.cmd.select('smasked', smasked)
		pymol.cmd.select('qpresent', qpresent)
		pymol.cmd.select('spresent', spresent)
		pymol.cmd.select('qaligned', qaligned)
		pymol.cmd.select('saligned', saligned)
		pymol.cmd.hide('everything')
		pymol.cmd.show('cartoon', '({}) or ({})'.format(qmasked, qpresent))
		pymol.cmd.show('cartoon', '({}) or ({})'.format(smasked, spresent))
		pymol.cmd.color('smudge', query)
		pymol.cmd.color('lightteal', subject)
		pymol.cmd.color('forest', 'qchain')
		pymol.cmd.color('deepteal', 'schain')
		pymol.cmd.color('palegreen', 'qpresent')
		pymol.cmd.color('palecyan', 'spresent')
		pymol.cmd.color('green', 'qaligned')
		pymol.cmd.color('cyan', 'saligned')
		pymol.cmd.center('qpresent or spresent')

	pymol.cmd.extend('spview', spview)
	spview()
	pymol.cmd.deselect()
	for fn in os.listdir(prefix):
		os.remove('{}/{}'.format(prefix, fn))
	os.rmdir(prefix)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile')
	parser.add_argument('name')

	args = parser.parse_args()

	main(args.name, args.infile)
