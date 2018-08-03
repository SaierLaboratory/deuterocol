#!/usr/bin/env python
from __future__ import print_function, division

import urllib
import sys, os, re
import argparse
import json

#TODO: move common classes out
import pdbtm_dbtool

class OPM(object):
	def __init__(self):
		pass

def get_database(prefix='.'):
	print('Fetching database...', file=sys.stderr)
	db = urllib.urlopen('http://opm.phar.umich.edu/OPM-2018-07-26.json')
	print('Saving database...', file=sys.stderr)
	with open('{}/opmall.json'.format(prefix), 'w') as f:
		f.write(db.read())

def build_database(fn, prefix):
	print('Unpacking database...', file=sys.stderr)

	db = {'atlas':[], 'citation':[], 'class':[], 'classification':[], 'family':[], 'history':[], 'membrane':[], 'protein':[], 'relatedproteins':[], 'relatedproteinstemp':[], 'relatedproteinstempnew':[], 'species':[], 'subunits':[], 'superfamily':[], 'type':[]}
	with open(fn) as f:
		n = 0
		for l in f:
			if not l.startswith('['): continue
			db[sorted(db)[n]] = json.loads(l)
			n += 1

	proteins = {}
	for protein in db['protein']:
		pdbid = protein['pdbid']
		proteins[pdbid] = pdbtm_dbtool.PDB()
		proteins[pdbid].label = pdbid

	related_proteins = {}
	parent_proteins = {}
	for relprot in db['relatedproteins']:
		rel_pdbid = relprot['pdbid']
		par_pdbid = relprot['protein_pdbid']

		if par_pdbid in proteins and rel_pdbid in proteins: pass
		elif par_pdbid in proteins:
			parent_proteins[rel_pdbid] = par_pdbid
			try: related_proteins[par_pdbid].append(rel_pdbid)
			except KeyError: related_proteins[par_pdbid] = [rel_pdbid]
			#proteins[rel_pdbid] = pdbtm_dbtool.PDB()
			#proteins[rel_pdbid].label = rel_pdbid
		#necessary because 1ym6 is misannotated as being a relative of 1okc
		else:
			parent_proteins[par_pdbid] = rel_pdbid
			try: related_proteins[rel_pdbid].append(par_pdbid)
			except KeyError: related_proteins[rel_pdbid] = [par_pdbid]
			#proteins[par_pdbid] = pdbtm_dbtool.PDB()
			#proteins[par_pdbid].label = par_pdbid

	#print(related_proteins['1ym6'])
	#for p in related_proteins: print(p, related_proteins[p])
	#for p in proteins: print(p)
	#exit()

	print('Writing assignments to {}/ASSIGNMENTS.TSV...'.format(prefix), file=sys.stderr)
	writeme = []
	for subunit in db['subunits']:
		pdbid = subunit['pdbid']
		letter = subunit['letter']
		#print(subunit['pdbid'], subunit['letter'], subunit['pdbid'] in proteins)
		if pdbid not in proteins:
			#necessary because 2iqo_? is an orphaned subunit
			if pdbid in parent_proteins: continue
			else: proteins[pdbid] = pdbtm_dbtool.PDB()
		chain = pdbtm_dbtool.Chain(label=letter, parent=proteins[pdbid])

		proteins[pdbid].chains.append(chain)
		segments = subunit['segment']
		spans = re.split(r'\),?\s*[0-9][0-9]?\s*\(?', segments)
		segstr = ''
		trouble = '2he6'
		for i in spans:
			s = i
			if pdbid == trouble: print('start:', s)
			s = re.sub(r'\n[A-Za-z]+', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\r?\n', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'^\s*[0-9]+\(\s*', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\s+', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'[0-9]+\(', '', s)
			if pdbid == trouble: print(s)
			s = re.sub(r'\)[^\-]*', '', s)
			if pdbid == trouble: print(s)
			if s.strip(): segstr += s + ','
			if pdbid == trouble: print('final:', s)
		#print(segments)
		#print(segstr[:-1])
		#print(sorted(db))
		writeme.append('{}_{}\t{}\n'.format(pdbid, letter, segstr[:-1]))
		if pdbid == trouble: print(writeme[-1])
		if pdbid in related_proteins:
			for relpdbid in related_proteins[pdbid]:
				writeme.append('{}_{}\t{}\n'.format(relpdbid, letter, segstr[:-1]))

	last100 = []
	with open('{}/ASSIGNMENTS.TSV'.format(prefix), 'w') as f:
		for s in sorted(writeme):
			if not last100: pass
			elif hash(s) in last100: continue

			f.write(s)
			last100.append(hash(s))
			if len(last100) > 100: last100.pop(0)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Manages OPM databases. Automatically fetches the OPM database if no options are specified. Run without any arguments, opm_dbtool will retrieve the OPM database, store i in opm, and unpack it.')

	parser.add_argument('-d', '--db', default='opmall.json', help='name of concatenated database file (default:opmall.json)')
	parser.add_argument('-b', '--build-db', action='store_true', help='(re)build database from an existing opmall.json file (available at http://opm.phar.umich.edu/OPM-2018-07-26.json)')
	parser.add_argument('directory', nargs='?', default='opm', help='directory to store database in')
	parser.add_argument('-f', '--force-refresh', action='store_true', help='force overwrite of existing database. Functionally equivalent to removing the old database and rerunning.')
	parser.add_argument('-x', '--extract-sequences', action='store_true', help='extract sequences from OPM database too')

	args = parser.parse_args()

	if args.build_db: build_database(args.db, args.directory)
	else:
		if not os.path.isdir(args.directory): os.mkdir(args.directory)
		if args.force_refresh or not os.path.isfile('{}/{}'.format(args.directory, args.db)): get_database(args.directory)
		build_database('{}/{}'.format(args.directory, args.db), args.directory)

	#if args.extract_sequences:

