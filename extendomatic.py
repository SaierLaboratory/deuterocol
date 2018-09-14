#!/usr/bin/env python

from __future__ import print_function, division, generators
import argparse, json
import numpy as np
import scipy.stats
import sys

#import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
#matplotlib.use('TkAgg')

def info(*things): print('[INFO]:', *things, file=sys.stderr)

def unpack_obj(obj):
	rmsd = obj['rmsd']
	length = obj['length']

	aligned = 0
	for dist in obj['distances']: aligned += 0 if dist is None else 1
	qpresent, spresent = obj['qpresent'], obj['spresent']
	qpresent, spresent = qpresent[0][1] - qpresent[0][0] + 1, spresent[0][1] - spresent[0][0] + 1

	covs = aligned/qpresent, aligned/spresent
	mincov = min(covs)
	maxcov = max(covs)

	distances = [np.nan if x is None else x for x in obj['distances']]
	#print(np.nanmax(distances))

	return rmsd, length, mincov, maxcov

class Dataset(object):
	def __init__(self, f, count=1000, mode=None, marg=None):
		self.names = []
		self.rmsds = []
		self.lengths = []
		self.mincovs = []
		self.maxcovs = []
		self.count = count

		self.kernel = None
		self.pdf = None

		if mode == 'onebest': 
			info('Using selection criterion ONEBEST')
			self.parse_best_data(f)
		elif mode == 'stretchone': 
			info('Using selection criterion STRETCHONE')
			self.parse_best_data(f)
			f.seek(0)
			self.stretch_data(f, n=marg)
		else: 
			info('Using selection criterion ANY')
			self.parse_data(f)

	def get_dict(self):
		obj = {}
		for name, rmsd, length, mincov, maxcov in zip(self.names, self.rmsds, self.lengths, self.mincovs, self.maxcovs):
			obj[name] = {'rmsd':rmsd, 'length':length, 'mincov':mincov, 'maxcov':maxcov}
		return obj

	def stretch_data(self, f, n=1):
		n = 1 if n is None else n

		names = self.names[:]

		seeds = {}
		for name in names:
			query, qchain, qhel, vs, subject, schain, shel = name.split('_')
			qpdbc = '{}_{}'.format(query, qchain)
			spdbc = '{}_{}'.format(subject, schain)
			qhel = [int(x) for x in qhel[1:].split('-')]
			shel = [int(x) for x in shel[1:].split('-')]
			seeds[(qpdbc, spdbc)] = [qhel, shel]


		names = []
		rmsds = []
		lengths = []
		mincovs = []
		maxcovs = []
		best = {}
		for k in seeds: best[k] = []

		f.seek(0)
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue
			
			sl = l.split('\t')
			query, qchain, qhel, vs, subject, schain, shel = sl[0].split('_')
			qpdbc = '{}_{}'.format(query, qchain)
			spdbc = '{}_{}'.format(subject, schain)
			qhel = [int(x) for x in qhel[1:].split('-')]
			shel = [int(x) for x in shel[1:].split('-')]
			try:
				seedqhel, seedshel = seeds[(qpdbc, spdbc)]
				if qhel[0] <= seedqhel[0] and shel[0] >= seedshel[0]: continue
				elif qhel[0] >= seedqhel[0] and shel[0] <= seedshel[0]: continue
				elif qhel[0] <= seedqhel[0] <= qhel[1]: continue
				elif qhel[0] <= seedqhel[1] <= qhel[1]: continue
				elif seedqhel[0] <= qhel[0] <= seedqhel[1]: continue
				elif seedqhel[0] <= qhel[1] <= seedqhel[1]: continue
				elif shel[0] <= seedshel[0] <= shel[1]: continue
				elif shel[0] <= seedshel[1] <= shel[1]: continue
				elif seedshel[0] <= shel[0] <= seedshel[1]: continue
				elif seedshel[0] <= shel[1] <= seedshel[1]: continue
			except KeyError: continue

			obj = json.loads(sl[1])
			rmsd, length, mincov, maxcov = unpack_obj(obj)
			best[(qpdbc, spdbc)].append((mincov, rmsd, name, length, maxcov))

		for k in best:
			if not best[k]: continue
			for extra in sorted(best[k])[::-1][:n]:
				mincov, rmsd, name, length, maxcov = extra
				names.append(name)
				rmsds.append(rmsd)
				lengths.append(length)
				mincovs.append(mincov)
				maxcovs.append(maxcov)

		self.names = np.hstack([self.names, names])
		self.rmsds = np.hstack([self.rmsds, rmsds])
		self.lengths = np.hstack([self.lengths, lengths])
		self.mincovs = np.hstack([self.mincovs, mincovs])
		self.maxcovs = np.hstack([self.maxcovs, maxcovs])

	def parse_best_data(self, f):
		n = 0

		bestkeys = []
		best = {}

		names = []
		rmsds = []
		lengths = []
		mincovs = []
		maxcovs = []

		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue

			sl = l.split('\t')
			name = sl[0]
			obj = json.loads(sl[1])
			qpdbc = '{0}_{1}'.format(*sl[0].split('_'))
			spdbc = '{4}_{5}'.format(*sl[0].split('_'))
			rmsd, length, mincov, maxcov = unpack_obj(obj)
			if rmsd == -1: continue

			try:
				i = best[(qpdbc, spdbc)]
				if mincov > mincovs[i]:
					names[i] = name
					rmsds[i] = rmsd
					lengths[i] = length
					mincovs[i] = mincov
					maxcovs[i] = maxcov
				elif mincov == mincovs[i] and rmsd < rmsds[i]:
					names[i] = name
					rmsds[i] = rmsd
					lengths[i] = length
					mincovs[i] = mincov
					maxcovs[i] = maxcov
				else: pass
			except KeyError:
				names.append(name)
				rmsds.append(rmsd)
				lengths.append(length)
				mincovs.append(mincov)
				maxcovs.append(maxcov)
				best[(qpdbc, spdbc)] = len(mincovs) - 1
				bestkeys.append((qpdbc, spdbc))

			if len(best) > 300: 
				best.pop(bestkeys.pop(0))
			#print('best', len(best), 'bestkeys', len(bestkeys), 'mincovs', len(mincovs))
			

			self.names = np.array(names)
			self.rmsds = np.array(rmsds)
			self.lengths = np.array(lengths)
			self.mincovs = np.array(mincovs)
			self.maxcovs = np.array(maxcovs)
			#self.rmsds.append(rmsd)
			#self.lengths.append(length)
			#self.mincovs.append(mincov)
			#self.maxcovs.append(maxcov)

			n += 1
			if len(mincovs) == self.count: break
		#print('bestkeys:', len(bestkeys))

	def parse_data(self, f):
		n = 0
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue

			sl = l.split('\t')
			self.names.append(sl[0])
			obj = json.loads(sl[1])
			rmsd, length, mincov, maxcov = unpack_obj(obj)
			if rmsd == -1: continue

			self.rmsds.append(rmsd)
			self.lengths.append(length)
			self.mincovs.append(mincov)
			self.maxcovs.append(maxcov)

			n += 1
			if n == self.count: break

		self.names = np.array(self.names)
		self.rmsds = np.array(self.rmsds)
		self.lengths = np.array(self.lengths)
		self.mincovs = np.array(self.mincovs)
		self.maxcovs = np.array(self.maxcovs)

	def evaluate(self, *args, **kwargs): return self.kernel.evaluate(*args, **kwargs)

	def gen_rmsd_mincov_kde(self, min_rmsd=0.0, max_rmsd=7.0, min_mincov=0.0, max_mincov=1.0, rmsd_resolution=100, mincov_resolution=100):
		rmsds, mincovs = np.mgrid[min_rmsd:max_rmsd:rmsd_resolution*1j, \
			min_mincov:max_mincov:mincov_resolution*1j]
		positions = np.vstack([rmsds.ravel(), mincovs.ravel()])
		values = np.vstack([self.rmsds, self.mincovs])
		self.kernel = scipy.stats.gaussian_kde(values)
		self.pdf = np.reshape(self.kernel(positions).T, rmsds.shape)
		
class Plot(object):
	def __init__(self, fig=None, canvas=None, ax=None):
		self.fig = Figure() if fig is None else fig
		self.canvas = FigureCanvas(self.fig) if fig is None else canvas
		self.ax = self.fig.add_subplot(111) if ax is None else ax
		
def plot_kde(kde, fig, ax, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100):
	X, Y = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([X.ravel(), Y.ravel()])
	Z = np.reshape(kde.evaluate(positions).T, X.shape)

	aspect = (rmsdlim[1] - rmsdlim[0]) / (mincovlim[1] - mincovlim[0]) / 1.414
	im = ax.imshow(np.rot90(Z), cmap='magma', extent=rmsdlim+mincovlim, aspect=aspect)
	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')

	fig.colorbar(im, ax=ax)

def plot_densities(positive, negative, unknown, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100, density_outfile='kde_plot.png', scatter_outfile='kde_scatter.png'):
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []
	plots = []
	for i in range(3): 
		ax.append(figure.add_subplot(2, 2, len(ax)+1))
		plots.append(Plot(fig=figure, canvas=canvas, ax=ax[-1]))
		#print(dir(ax[-1]))

	plot_kde(positive, fig=figure, ax=ax[0], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[0].set_title('Positive (n={})'.format(len(positive.mincovs)))
	plot_kde(negative, fig=figure, ax=ax[1], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[1].set_title('Negative (n={})'.format(len(negative.mincovs)))
	plot_kde(unknown, fig=figure, ax=ax[2], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[2].set_title('Unknown (n={})'.format(len(unknown.mincovs)))

	figure.savefig(density_outfile, dpi=600)
	plot_rmsd_cov(positive, fig=figure, ax=ax[0])
	plot_rmsd_cov(negative, fig=figure, ax=ax[1])
	plot_rmsd_cov(unknown, fig=figure, ax=ax[2])
	figure.savefig(scatter_outfile, dpi=600)

def plot_rmsd_cov(kde, fig, ax):
	ax.plot(kde.rmsds, kde.mincovs, c=(1., 1., 1., 0.1), marker='.', linewidth=0)

def main(positivefn, negativefn, unknownfn, count=1000, density_outfile='density_plot.png', scatter_outfile='scatter_plot.png', stretch=1):
	with open(positivefn) as f: positive = Dataset(f, count=count, mode='stretchone', marg=stretch)
	positive.gen_rmsd_mincov_kde()
	with open(negativefn) as f: negative = Dataset(f, count=count, mode='stretchone', marg=stretch)
	negative.gen_rmsd_mincov_kde()
	with open(unknownfn) as f: unknown = Dataset(f, count=count, mode='stretchone', marg=stretch)
	unknown.gen_rmsd_mincov_kde()

	print('#{}, {}, {}. stretch={}, n={}'.format(positivefn, negativefn, unknownfn, stretch, count))

	min_rmsd, max_rmsd = 0.0, 6.0
	min_mincov, max_mincov = 0.4, 1.0

	unkpoints = np.vstack([unknown.rmsds, unknown.mincovs])
	unkposteriors = positive.evaluate(unkpoints) * len(positive.rmsds) / (positive.evaluate(unkpoints) * len(positive.rmsds) + negative.evaluate(unkpoints) * len(negative.rmsds))
	print('#Alignment\tPosterior\tRMSD\tCoverage')
	for post, name, rmsd, mincov in sorted(zip(unkposteriors, unknown.names, unknown.rmsds, unknown.mincovs))[::-1][:10]:
		print('{}\t{:0.02e}\t{:0.1f}\t{:0.0%}'.format(name, post, rmsd, mincov))
	plot_densities(positive, negative, unknown, rmsdlim=(min_rmsd, max_rmsd), mincovlim=(0.0, max_mincov), resolution=100, density_outfile=density_outfile, scatter_outfile=scatter_outfile)

	#print('stats: ({} <= RMSD <= {}, {} <= minCov <= {}, n={})'.format(min_rmsd, max_rmsd, min_mincov, max_mincov, count*3))
	#print('positive:', positive.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))
	#print('negative:', negative.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))
	#print('unknown:', unknown.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('positive')
	parser.add_argument('negative')
	parser.add_argument('unknown')
	parser.add_argument('-n', type=int, default=1000)
	parser.add_argument('-o', nargs=2, default=['density_plot.png', 'scatter_plot.png'])
	parser.add_argument('-s', default=1, type=int)

	args = parser.parse_args()

	main(args.positive, args.negative, args.unknown, count=args.n, density_outfile=args.o[0], stretch=args.s, scatter_outfile=args.o[1])
