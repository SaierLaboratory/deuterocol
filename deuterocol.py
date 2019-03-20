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
		elif params.prop == 'tmdata': self.tmdatadir = params.val[0]
		#2 vals required

	def get_children(self, tclist, level=4):
		alltcdb = TCIDCollection()
		with open('{}/tcmap.tsv'.format(self.tmdatadir)) as f:
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

	def run(self):
		if self.children is not None:
			self.fams1 = [str(fam) for fam in self.get_children(self.fams1, level=self.children)]
			self.fams2 = [str(fam) for fam in self.get_children(self.fams2, level=self.children)]

		if not os.path.isdir(self.outdir): os.mkdir(self.outdir)
		deut1 = Deuterocol1.Deuterocol1(tmdatadir=self.tmdatadir, outdir='{}/deuterocol1'.format(self.outdir), 
			invert=self.generate_negative, inclusive=self.generate_negative, 
		)
		deut1.run(*(self.fams1 + self.fams2), force=self.force)

		deut2 = Deuterocol2.Deuterocol2(tmdatadir=self.tmdatadir, d1dir='{}/deuterocol1'.format(self.outdir),
			outdir=self.outdir, bundle=self.bundle, allow_internal=self.allow_internal
		)
		deut2.run(self.fams1, self.fams2)

		tma = tmalign.TMalign(d2dir=self.outdir, skip_cut=self.skip_cut)
		tma.run()
		#TODO: min_tms, max_tms

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('--config', help='Load config file')
	parser.add_argument('outdir', nargs='?')
	args = parser.parse_args()

	deuterocol = Deuterocol() if not args.outdir else Deuterocol(outdir=args.outdir)
	if args.config: 
		with open(args.config) as f: deuterocol.parse_config(f.read())
	deuterocol.run()
