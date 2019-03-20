#!/usr/bin/env python
from __future__ import print_function, division

import urllib
import sys, os, re
import argparse
import json
import math
import csv
import subprocess

import Bio.PDB
try: CODE = Bio.PDB.protein_letters_3to1
except AttributeError: CODE = Bio.PDB.to_one_letter_code

#TODO: move common classes out
import pdbtm_dbtool

class OPM(object):
	def __init__(self):
		pass


def _scroll_json(url):
	''' Collects all pages for a multipage JSON '''
	out = []
	f = urllib.urlopen(url)
	raw = f.read()
	obj = json.loads(raw)
	total_objects = obj['total_objects']
	page_end = obj['page_end']
	page_num = obj['page_num']
	page_size = obj['page_size']

	out += obj['objects']

	for i in range(page_num+1, int(math.ceil(total_objects/page_size))):
		f = urllib.urlopen(url + '?pageNum={}'.format(i))
		raw = f.read()
		obj = json.loads(raw)
		out += obj['objects']
	return json.dumps(out)
def get_database(prefix='.'):
	''' (Deprecated) Fetched the OPM database dump '''
	print('Fetching database...', file=sys.stderr)
	with open('{}/opmall.json'.format(prefix), 'w') as f: pass

	for table in ('types', 'classtypes', 'superfamilies', 'families', 'primary_structures', 'structure_subunits'):
		data = _scroll_json('https://lomize-group-opm.herokuapp.com/{}'.format(table))
		
		#db = urllib.urlopen('https://lomize-group-opm.herokuapp.com/{}'.format(table))
		with open('{}/opmall.json'.format(prefix), 'a') as f: 
			f.write('// {}\n'.format(table))
			f.write(data)
			f.write('\n')
	#print('Saving database...', file=sys.stderr)
	exit()

def build_database(fn, prefix):
	''' Unpacked the OPM database and converted it into the ASSIGNMENTS.TSV '''
	print('Unpacking database...', file=sys.stderr)

	#db = {'atlas':[], 'citation':[], 'class':[], 'classification':[], 'family':[], 'history':[], 'membrane':[], 'protein':[], 'relatedproteins':[], 'relatedproteinstemp':[], 'relatedproteinstempnew':[], 'species':[], 'subunits':[], 'superfamily':[], 'type':[]}
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

def correct_spans(segments, pdbid=None):
	''' Corrected some of the easier formatting issues. TODO: handle re-entrants properly '''

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
	return segstr[:-1]

	writeme.append('{}_{}\t{}\n'.format(pdbid, letter, segstr[:-1]))

def process_opm_csv(outdir):
	''' Parsed OPM CSVs '''
	if not os.path.isdir(outdir): os.mkdir(outdir)

	try: urlopen = urllib.urlopen
	except AttributeError: urlopen = urllib.request.urlopen

	if 0:
		inf = urlopen('https://lomize-group-opm.herokuapp.com/structure_subunits?fileFormat=csv')
		with open('{}/opm_subunits.csv'.format(outdir), 'w') as outf:
			outf.write(inf.read())
		inf.close()

	spans = {}
	with open('{}/opm_subunits.csv'.format(outdir)) as f:
		csv_reader = csv.reader(f)
		for row in csv_reader:
			pdbid = re.findall('[0-9a-zA-Z]+', row[5])[0].lower()
			chain = row[6]
			segments = row[11]
			print(correct_spans(segments, pdbid=pdbid))
	

def process_opm_scroll_primaries(outdir, force=False):
	''' Processed primary structure records from OPM, linking secondary representations to their primaries 

	TODO: Put everything in a class and break up the easily decoupled steps into separate methods
	'''
	if not os.path.isdir(outdir): os.mkdir(outdir)

	try: urlopen = urllib.urlopen
	except AttributeError: urlopen = urllib.request.urlopen

	if force or not os.path.isfile('{}/primary_structures.json'.format(outdir)):
		primary_structures = _scroll_json('https://lomize-group-opm.herokuapp.com/primary_structures')
		with open('{}/primary_structures.json'.format(outdir), 'w') as f:
			f.write(primary_structures)
	else: print('primary_structures.json already exists, skipping download')

	idlist = []
	with open('{}/primary_structures.json'.format(outdir)) as f:
		struclist = json.loads(f.read())
		for struc in struclist: idlist.append(struc['id'])
	
	primary_structures = {}
	if force or not os.path.isfile('{}/indiv_structures.tsv'.format(outdir)):
		with open('{}/indiv_structures.tsv'.format(outdir), 'w') as f:
			for n, opmid in enumerate(idlist):
				url = urlopen('https://lomize-group-opm.herokuapp.com/primary_structures/{}'.format(opmid))
				#print(json.loads(url.read()))
				raw = url.read()
				obj = json.loads(raw)
				primary_structures[obj['pdbid']] = obj
				f.write('{}\t{}\n'.format(obj['pdbid'], raw))

				if not (n % 100): 
					f.flush()
					print('Committed {}/{} individual structure records to disk'.format(n, len(idlist)))
	else: print('indiv_structures.json already exists, skipping download')

	secondary_representations = {}
	with open('{}/indiv_structures.tsv'.format(outdir)) as f:
		for l in f:
			if not l.strip(): continue
			elif l.lstrip().startswith('#'): continue

			pdbid, raw = l.split('\t')
			obj = json.loads(raw)
			primary_structures[pdbid] = obj
			if obj['secondary_representations']:
				secondary_representations[pdbid] = obj['secondary_representations']

	for primaryid in secondary_representations:
		for secondaryobj in secondary_representations[primaryid]:
			secondaryid = secondaryobj['pdbid']
			if secondaryid in primary_structures: continue

			else: 
				primary_structures[secondaryid] = primary_structures[primaryid].copy()
				primary_structures[secondaryid]['pdbid'] = secondaryid
				#print(primary_structures[primaryid]['pdbid'], primary_structures[secondaryid]['pdbid'])

	if force or not os.path.isfile('{}/ASSIGNMENTS.TSV'.format(outdir)):
		with open('{}/ASSIGNMENTS.TSV'.format(outdir), 'w') as f:
			for pdbid in sorted(primary_structures):
				for subunit in primary_structures[pdbid]['subunits']:
					chain = subunit['protein_letter']
					segment = subunit['segment']

					spans = correct_spans(segment, pdbid=pdbid)
					f.write('{}_{}\t{}\n'.format(pdbid.lower(), chain, spans))
	else: print('ASSIGNMENTS.TSV already exists, skipping rebuild')

	print('Mapping secondaries back to primaries')
	sec2prim = {}
	for primaryid in secondary_representations:
		for secobj in secondary_representations[primaryid]:
			sec2prim[secobj['pdbid']] = primaryid


	if force or not os.path.isfile('{}/all_pdbs.tar.gz'.format(outdir)):
		subprocess.call(['wget', '-O', '{}/all_pdbs.tar.gz'.format(outdir), 'https://storage.googleapis.com/opm-assets/pdb/tar_files/all_pdbs.tar.gz'])
	else: print('Found all_pdbs.tar.gz, skipping redownload')

	if force or not os.path.isdir('{}/pdb'.format(outdir)): 
		if not os.path.isdir('{}/pdb'.format(outdir)): os.mkdir('{}/pdb'.format(outdir))

		subprocess.call(['tar', 'xzf', '{}/all_pdbs.tar.gz'.format(outdir), '-C', outdir])
		#f_remote = urlopen('https://storage.googleapis.com/opm-assets/pdb/tar_files/all_pdbs.tar.gz', 'rb')
		#with open('{}/all_pdbs.tar.gz'.format(outdir), 'wb') as f_local:
		#	f_local.write(f_remote.read())
	else: print('Found pdb subdirectory, skipping reextraction')


	if not os.path.isdir('{}/sequences'.format(outdir)): os.mkdir('{}/sequences'.format(outdir))


	sequences = {}
	if not os.path.isdir('{}/sequences'.format(outdir)): os.mkdir('{}/sequences'.format(outdir))

	print('Counting relevant PDB subunits')
	relevant_pdbc = set()
	for pdbid in sorted(primary_structures):
		for subunit in primary_structures[pdbid]['subunits']:
			chain = subunit['protein_letter']
			pdbc = '{}_{}'.format(pdbid.upper(), chain)
			relevant_pdbc.add(pdbc)

	print('Extracting {} sequences. This may take some time'.format(len(primary_structures)))
	for n, pdbid in enumerate(sorted(primary_structures)):
		if not n % 500: print('Extracted {} sequences so far'.format(n))
		fn = '{}/pdb/{}.pdb'.format(outdir, pdbid.lower())
		if not os.path.isfile(fn):
			if pdbid in sec2prim: fn = '{}/pdb/{}.pdb'.format(outdir, sec2prim[pdbid].lower())
			else: 
				subprocess.call(['wget', '-O', fn + '.gz', 'https://files.rcsb.org/download/{}.pdb.gz'.format(pdbid)])
				subprocess.call(['gunzip', fn + '.gz'])

		done_seq = True
		for subunit in primary_structures[pdbid]['subunits']:
			chain = subunit['protein_letter']

			if force or not os.path.isfile('{}/sequences/{}_{}.fa'.format(outdir, pdbid, chain)):
				done_seq = False
				break
		if done_seq: continue

		if os.path.isfile(fn):
			with open(fn) as f:
				pdbseq = extract_pdb_sequences(f)
			for chain in pdbseq:
				pdbc = '{}_{}'.format(pdbid.upper(), chain)
				if pdbc not in relevant_pdbc: continue

				outfn = '{}/sequences/{}_{}.fa'.format(outdir, pdbid.upper(), chain)
				if force or not os.path.isfile(outfn):
					with open(outfn, 'w') as f:
						f.write('>{}_{}\n{}'.format(pdbid.upper(), chain, pdbseq[chain]))


	tcdb2pdb = {}

	print('Evaluating sequence quality')
	skipme = set()
	if os.path.isfile('{}/bad_seqlist'.format(outdir)): 
		with open('{}/bad_seqlist'.format(outdir)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.lstrip().startswith('#'): continue
				skipme.add(l.strip())

	skipf = open('{}/bad_seqlist'.format(outdir), 'w')
	print('BLASTing sequences for {} structures against TCDB. This could take a while.'.format(len(primary_structures)))
	for n, pdbid in enumerate(sorted(primary_structures)):
		if not n % 500: print('BLASTed sequences for {} structures so far'.format(n))
		for subunit in primary_structures[pdbid]['subunits']:
			chain = subunit['protein_letter']
			pdbc = '{}_{}'.format(pdbid.upper(), chain)

			fn = '{}/sequences/{}.fa'.format(outdir, pdbc)

			if skipme and (pdbc in skipme): continue
			else:
				if os.path.isfile(fn):
					with open(fn) as f:
						if is_low_quality_sequence(f.read()): 
							skipf.write(pdbc + '\n')
							continue
					tcid = quickblast(fn, expect=1e-5, ident=0.95)
					if tcid: 
						try: tcdb2pdb[tcid].append(pdbc)
						except KeyError: tcdb2pdb[tcid] = [pdbc]

	if force or not os.path.isfile('{}/tcmap.tsv'.format(outdir)):
		with open('{}/tcmap.tsv'.format(outdir), 'w') as f:
			line = ''
			for tcid in sorted(tcdb2pdb):
				line += tcid + '\t'
				for pdb in sorted(tcdb2pdb[tcid]):
					line += pdb + ','
			f.write(line[:-1] + '\n')
	else: print('OPM-specific tcmap found, skipping remapping')
				

def extract_pdb_sequences(f):
	''' Grabs sequences from the ATOM records in a PDB file '''
	sequences = {}
	for l in f:
		if not l.strip(): continue
		elif not l.startswith('ATOM'): continue
		if l[13:15] != 'CA': continue
		try: resn = CODE[l[17:20]]
		except KeyError: continue
		#resi = l[22:26]
		chain = l[21]
		if not chain.strip(): chain = 'A'

		try: sequences[chain] += resn
		except KeyError: sequences[chain] = resn
	return sequences
		
def is_low_quality_sequence(fasta):
	''' Heuristic for a sequence being too incomplete to even BLAST '''
	seq = re.sub('[^A-Za-z*]', '', fasta[fasta.find('\n'):])
	if seq.count('X') > (len(seq)*.75): return True
	elif len(seq) < 20: return True
	return False


def quickblast(fn, expect=1e-5, ident=0.95):
	''' BLAST a sequence against TCDB. Requires that the TCDB BLAST database be somewhere in $BLASTDB '''
	out = subprocess.check_output(['blastp', '-db', 'tcdb', '-outfmt', '6', '-evalue', str(expect), '-max_target_seqs', '1', '-query', fn])
	if not out.strip(): return None

	row = out.split('\t')
	e = float(row[10])
	i = row[2]
	acc = row[1]


	if i >= ident: return acc
	elif e <= expect: return acc
	else: return None


def fetch_seq(pdbc):
	''' Grabs a sequence from PDBAA '''
	p = subprocess.Popen(['blastdbcmd', '-db', 'pdbaa', '-entry', pdbc, '-target_only'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = p.communicate()

	return out

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='Manages OPM databases. Automatically fetches the OPM database if no options are specified. Run without any arguments, opm_dbtool will retrieve the OPM database, store i in opm, and unpack it.')

	parser.add_argument('directory', nargs='?', default='opm', help='directory to store database in')
	parser.add_argument('-f', '--force-refresh', action='store_true', help='force overwrite of existing database. Functionally equivalent to removing the old database and rerunning.')

	args = parser.parse_args()

	process_opm_scroll_primaries(outdir=args.directory, force=args.force_refresh)
