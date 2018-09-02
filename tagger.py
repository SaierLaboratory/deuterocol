#!/usr/bin/env python2
from __future__ import print_function, division, generators

import xml.etree.ElementTree as ET
import argparse
import pdbtm_dbtool

def tag(f):
	pdb = pdbtm_dbtool.PDB.parse_xml(f)

	for chain in pdb:
		if not chain.count(): continue
		#print(chain.seq)
		for region in chain: 
			if region.topo in pdbtm_dbtool.TRANSMEMBRANE: 
				print(region.seq)

	return 0

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('infile', nargs='+')


	args = parser.parse_args()

	for fn in args.infile:
		with open(fn) as f:
			print(tag(f))

