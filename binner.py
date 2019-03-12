#!/usr/bin/env python2
from __future__ import print_function, division

import scipy.stats, json, argparse
import numpy as np

import matplotlib
#matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter


def unpack_spans(l):
	out = []
	for span in l:
		if span[0] is None: out += [None] * span[1]
		else: out += list(range(span[0], span[1]+1))
	return out
def extract_info(obj):
	qpresent = 0
	for span in obj['qpresent']: qpresent += span[1] - span[0] + 1
	spresent = 0
	for span in obj['spresent']: spresent += span[1] - span[0] + 1

	qaligned = 0
	#for span in obj['qaligned']:
	#	if span[0] is None: continue
	#	else: qaligned += span[1] - span[0]
	#print(obj['qaligned'])
	#print(unpack_spans(obj['qaligned']))
	#print(len(unpack_spans(obj['qaligned'])))
	#print(len(unpack_spans(obj['saligned'])))
	#print(len(obj['distances']))
	saligned = 0
	#for span in obj['saligned']:
	#	if span[0] is None: continue
	#	else: saligned += span[1] - span[0]
	aligned = 0
	for dist in obj['distances']:
		if dist is None: continue
		else: aligned += 1

	#return obj['rmsd'], min(qaligned/qpresent, saligned/spresent)
	return obj['rmsd'], aligned/max(qpresent, spresent)

def get_bin(bins, data):
	for i, x in enumerate(bins):
		if data <= x: return i
	return i

def plot_kernel(kernel, lim=None, values=None):
	lim = [[0., 10.], [0., 1.]] if lim is None else lim
	#X, Y = np.mgrid[lim[0][0]:lim[0][1]:0.1, lim[1][0]:lim[1][1]:0.01]
	X, Y = np.mgrid[lim[0][0]:lim[0][1]:0.05, lim[1][0]:lim[1][1]:0.005]
	positions = np.vstack([X.ravel(), Y.ravel()])
	Z = np.reshape(kernel(positions).T, X.shape)
	
	fig, ax = plt.subplots()
	plt.tight_layout(True)
	#fig.set_figwidth(10)
	#fig.set_figheight(10)
	#ax.imshow(np.rot90(Z), cmap=cm.gist_earth_r, extent=lim[0]+lim[1], aspect=10)
	ax.imshow(np.rot90(Z), cmap=cm.magma, extent=lim[0]+lim[1], aspect=7)
	#ax.plot(values[0], values[1], '.', markersize=1)
	ax.set_xlim(lim[0])
	ax.set_ylim(lim[1])
	plt.show()
	
def plot(x, y, z):
	X, Y = np.meshgrid(x, y)
	plt.rc('text', usetex=True)
	fig = plt.figure()
	#ax = fig.add_subplot(111, projection='3d')
	#ax.plot_
	ax = fig.gca(projection='3d')
	#ax.contour(X, Y, z, extend3d=True, cmap=cm.coolwarm)
	surf = ax.plot_surface(X, Y, z/np.sum(z), antialiased=False, cmap=cm.magma)
	#ax.set_zscale('log')
	ax.set_zlim(0.0, np.max(z/np.sum(z)))
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.03f'))
	fig.colorbar(surf, shrink=0.5, aspect=5)


	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')
	ax.set_zlabel('f')

	plt.savefig('qr.png')
	plt.show()

def main(filelist, max_rmsd=4.0, min_cov=0.8):
	bins_rmsd = np.arange(0, max_rmsd, 0.1)
	bins_cov = np.arange(min_cov, 1.0, 0.01)

	freq_rmsd = np.zeros(len(bins_rmsd))
	freq_cov = np.zeros(len(bins_cov))
	freq_2d = np.zeros([len(freq_cov), len(freq_rmsd)])

	rmsds = []
	covs = []
	n = 0
	for fn in filelist:
		with open(fn) as f:
			for l in f: 
				sl = l.split('\t')
				name = sl[0]
				obj = json.loads(sl[1])
				rmsd, cov = extract_info(obj)

				#if rmsd > max_rmsd: continue
				if rmsd == -1: continue
				#if cov < min_cov: continue
				if cov == -1: continue

				rmsdbin = get_bin(bins_rmsd, rmsd)
				covbin = get_bin(bins_cov, cov)
				freq_rmsd[rmsdbin] += 1
				freq_cov[covbin] += 1
				freq_2d[covbin,rmsdbin] += 1
				

				rmsds.append(rmsd)
				covs.append(cov)
				#print(name, rmsd, cov)
	#print(zip(bins_rmsd, freq_rmsd))
	#print(zip(bins_cov, freq_cov))
	#print(freq_2d)
	values = np.vstack([rmsds, covs])
	kernel = scipy.stats.gaussian_kde(values)
	plot_kernel(lim=[[0., max_rmsd], [min_cov, 1.]], kernel=kernel, values=values)

	#plot(bins_rmsd, bins_cov, freq_2d)


if __name__ == '__main__':
	parser = argparse.ArgumentParser()

	parser.add_argument('infile', nargs='+', help='superpositions to do stats for')
	parser.add_argument('--max-rmsd', type=float, default=7.5)
	parser.add_argument('--min-cov', type=float, default=0.0)

	args = parser.parse_args()

	main(args.infile, max_rmsd=args.max_rmsd, min_cov=args.min_cov)
