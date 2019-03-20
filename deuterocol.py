#!/usr/bin/env python
#This script automatically runs Deuterocol1.py, Deuterocol2.py, and tmalign.py together, thus raising the question of why they ever needed to be separate to begin with

import argparse
import shlex
try: 
	import urllib.request
	urlopen = urllib.request.urlopen
except ImportError:
	import urllib
	urlopen = urllib.urlopen
import re
import os

import Deuterocol1
import Deuterocol2
import tmalign
import json2tsv

from kevtools_common.types import TCID
from kevtools_common.types import TCIDCollection


class Deuterocol(object):
	fams1 = []
	fams2 = []
	generate_negative = False
	outdir = 'deuterocol_out'
	bundle = 0
	children = None
	tmdata = 'tmdata'
	force = False
	allow_internal = False
	skip_cut = False
	min_tms = 3

	def __init__(self, outdir='deuterocol_out'):
		self.outdir = outdir

	def parse_config(self, raw):
		for l in raw.split('\n'):
			if not l.strip(): continue
			elif l.lstrip().startswith('#'): continue
			else:
				sl = shlex.split(l)
				if sl[0].lower() == 'set': self._set_prop(args=sl[1:])

	def _set_prop(self, args):
		parser = argparse.ArgumentParser()
		parser.add_argument('prop')
		parser.add_argument('val', nargs='*')
		params = parser.parse_args(args)

		#0 vals allowed:
		if params.prop == 'fams1': self.fams1 = params.val
		elif params.prop == 'fams2': 
			if params.val[0] == 'negative':
				self.fams2 = []
				self.generate_negative = True
			else: self.fams2 = params.val
		#TODO: allow setting these boolean props to negative/false
		elif params.prop == 'force': self.force = True
		elif params.prop == 'allow_internal': self.allow_internal = True
		elif params.prop == 'skip_cut': self.skip_cut = True
		#1 val required:
		elif len(params.val) < 1: raise ArgumentError('Not enough values supplied for positional argument val')
		elif params.prop == 'outdir': self.outdir = params.val[0]
		elif params.prop == 'bundle': self.bundle = int(params.val[0])
		elif params.prop == 'children': self.children = int(params.val[0])
		elif params.prop == 'tmdata': self.tmdata = params.val[0]
		elif params.prop == 'min_tms': self.min_tms = int(params.val[0])
		#2 vals required

	def get_children(self, tclist, level=4):
		alltcdb = TCIDCollection()
		with open('{}/tcmap.tsv'.format(self.tmdata)) as f:
			for l in f:
				if not l.strip(): continue
				elif l.lstrip().startswith('#'): continue
				alltcdb.add_tcid(l.split('\t')[0])
		#url = 'http://www.tcdb.org/public/tcdb'
		#f = urlopen(url)
		#raw = f.read()
		#f.close()
		#for l in raw.split('\n'):
		#	if l.startswith('>'): 
		#		s = l[1:].strip()
		#		s = re.sub(' \|', '\|', s).split()[0]
		#		splits = s.split('|')
		#		tc = splits[-1]
		#		acc = splits[-2]
		#		
		#		alltcdb.add_tcid('{}-{}'.format(tc, acc))

		parents = TCIDCollection()
		for fam in tclist:
			parents.add_tcid(fam)
		children = TCIDCollection()
		for tcid in alltcdb:
			for fam in parents:
				if tcid in fam:  
					children.add_tcid(tcid[:level])
					break
		return sorted(children)

	def write_configuration(self):
		with open('{}/deuterocol.cfg'.format(self.outdir), 'w') as f:
			s = 'set outdir {}\n'.format(self.outdir)
			f.write(s)

			s = 'set fams1'
			for fam in self.fams1: s += ' {}'.format(fam)
			s += '\n'
			f.write(s)

			s = 'set fams2'
			for fam in self.fams2: s += ' {}'.format(fam)
			s += '\n'
			f.write(s)

			s = 'set children {}\n'.format(self.children)
			f.write(s)

			s = 'set bundle {}\n'.format(self.bundle)
			f.write(s)

			s = 'set tmdata {}\n'.format(self.tmdata)
			f.write(s)

			s = 'set min_tms {}\n'.format(self.min_tms)
			f.write(s)

			if self.force:
				s = 'set force\n'
				f.write(s)
			if self.allow_internal:
				s = 'set allow_internal\n'
				f.write(s)
			if self.skip_cut:
				s = 'set skip_cut\n'
				f.write(s)
	def run(self):
		if self.children is not None:
			self.fams1 = [str(fam) for fam in self.get_children(self.fams1, level=self.children)]
			self.fams2 = [str(fam) for fam in self.get_children(self.fams2, level=self.children)]

		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)

		self.write_configuration()

		Deuterocol1.info('Running Deuterocol 1...')
		deut1 = Deuterocol1.Deuterocol1(tmdatadir=self.tmdata, outdir='{}/deuterocol1'.format(self.outdir), 
			invert=self.generate_negative, inclusive=self.generate_negative, 
		)
		deut1.run(*(self.fams1 + self.fams2), force=self.force)

		Deuterocol1.info('Running Deuterocol 2...')
		deut2 = Deuterocol2.Deuterocol2(tmdatadir=self.tmdata, d1dir='{}/deuterocol1'.format(self.outdir),
			outdir=self.outdir, bundle=self.bundle, allow_internal=self.allow_internal
		)
		deut2.run(self.fams1, self.fams2)

		Deuterocol1.info('Running TMalign...')
		tma = tmalign.TMalign(d2dir=self.outdir, skip_cut=self.skip_cut, min_tms=self.min_tms)
		tma.run()
		#TODO: min_tms, max_tms

		Deuterocol1.info('Generating tables...')
		self.generate_tables()

	def generate_tables(self):
		if not os.path.isdir('{}/tables'.format(self.outdir)): os.mkdir('{}/tables'.format(self.outdir))

		json2tsv.main(self.outdir, '{}/tables/all.tsv'.format(self.outdir))

		for subdir in os.listdir(self.outdir):
			path = '{}/{}'.format(self.outdir, subdir)
			if '_vs_' in subdir:
				if os.path.isfile(path + '/tmalignments/sp_all.tsv'):
					infn = path + '/tmalignments/sp_all.tsv'
					outfn = '{}/tables/{}.tsv'.format(self.outdir, subdir)
					json2tsv.main(infn, outfn)

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help='Load config file')
	parser.add_argument('--force', action='store_true', help='Overwrite files with wild abandon')
	parser.add_argument('--allow-internal', action='store_true', help='Allow within-TCID comparisons')
	parser.add_argument('--skip-cut', action='store_true', help='Skip cutting step. Only use this if all PDBs have been cut.')
	parser.add_argument('--bundle', type=int, default=None, help='Bundle size')
	parser.add_argument('--children', type=int, default=None, help='What TC-ID level to split families to')
	if 'TMDATA' in os.environ: parser.add_argument('--tmdata', default=os.environ['TMDATA'], help='Where tmdata is (default: $TMDATA == {})'.format(os.environ['TMDATA']))
	else: parser.add_argument('--tmdata', default=None, help='Where tmdata is (default: tmdata)')
	parser.add_argument('--min-tms', type=int, default=None, help='Minimum TMSs required to process an alignment')
	parser.add_argument('--fams1', nargs='+', help='First set of families')
	parser.add_argument('--fams2', nargs='+', help='Second set of families. Specify "negative" to generate a negative control (not reimplemented yet)')

	parser.add_argument('outdir', nargs='?')
	args = parser.parse_args()

	deuterocol = Deuterocol() if not args.outdir else Deuterocol(outdir=args.outdir)
	if args.config: 
		with open(args.config) as f: deuterocol.parse_config(f.read())

	if args.force: deuterocol.force = True
	if args.allow_internal: deuterocol.allow_internal = True
	if args.skip_cut: deuterocol.skip_cut = True
	if args.bundle: deuterocol.bundle = args.bundle
	if args.children: deuterocol.children = args.children
	#TODO: reconcile with cfg
	if args.tmdata: deuterocol.tmdata = args.tmdata
	if args.min_tms: deuterocol.min_tms = args.min_tms
	if args.fams1: deuterocol.fams1 = args.fams1
	if args.fams2: deuterocol.fams2 = args.fams2
	deuterocol.run()

