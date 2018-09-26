#!/usr/bin/env python

from __future__ import print_function, division, generators
from mpl_toolkits.mplot3d import Axes3D
import argparse, json
import numpy as np
import scipy.stats
import sys

PRIOR = 0.5
BEAUTYFACTOR = 8

#import matplotlib
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
#matplotlib.use('TkAgg')

def info(*things): print('[INFO]:', *things, file=sys.stderr)

def unpack_obj(obj):
	rmsd = obj['rmsd']
	length = obj['length']

	oldaligned = 0
	for dist in obj['distances']: oldaligned += 0 if dist is None else 1

	qpresent, spresent = obj['qpresent'], obj['spresent']

	qaligned = []
	saligned = []
	for span in obj['qaligned']:
		if span[0] is None: qaligned += [None] * span[1]
		else: qaligned += list(range(span[0], span[1]+1))
	for span in obj['saligned']:
		if span[0] is None: saligned += [None] * span[1]
		else: saligned += list(range(span[0], span[1]+1))

	aligned = 0
	qp = 0
	sp = 0
	for q, dist, s in zip(qaligned, obj['distances'], saligned):
		#print(qpresent, q, dist, s, spresent)
		if q is None: continue
		elif (q < qpresent[0][0]) or (q > qpresent[0][1]): 
			rmsd = -1
			continue
		else: qp += 1
		if dist is None: continue
		if s is None: continue
		elif (s < spresent[0][0]) or (s > spresent[0][1]): 
			rmsd = -1
			continue
		else: sp += 1
		aligned += 1
	#if oldaligned != aligned: print(oldaligned, aligned)

	qpresent, spresent = qpresent[0][1] - qpresent[0][0] + 1, spresent[0][1] - spresent[0][0] + 1

	#old calc
	#minpresent = min(qpresent, spresent)
	minpresent = min(qp, sp)

	#covs = aligned/qpresent, aligned/spresent
	if (qp == 0) or (sp == 0):
		covs = 0, 0
	else: covs = aligned/qp, aligned/sp
	mincov = min(covs)
	maxcov = max(covs)

	#distances = [np.nan if x is None else x for x in obj['distances']]
	#print(np.nanmax(distances))

	return rmsd, length, mincov, maxcov, minpresent

class Dataset(object):
	def __init__(self, f, count=1000, mode=None, marg=None, min_present=50):
		self.names = []
		self.rmsds = []
		self.lengths = []
		self.mincovs = []
		self.maxcovs = []
		self.count = count

		self.min_present = min_present
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

			try: obj = json.loads(sl[1])
			except ValueError: continue
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
			rmsd, length, mincov, maxcov, minpresent = unpack_obj(obj)
			if minpresent < self.min_present: continue
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
		names = []
		rmsds = []
		lengths = []
		mincovs = []
		maxcovs = []
		for l in f:
			if not l.strip(): continue
			elif l.startswith('#'): continue

			sl = l.split('\t')
			obj = json.loads(sl[1])
			rmsd, length, mincov, maxcov, minpresent = unpack_obj(obj)
			if minpresent < self.min_present: continue
			if rmsd == -1: continue
			if not mincov: continue

			names.append(sl[0])
			rmsds.append(rmsd)
			lengths.append(length)
			mincovs.append(mincov)
			maxcovs.append(maxcov)

			n += 1
			#if n == self.count: break

		self.names = []
		self.rmsds = []
		self.lengths = []
		self.mincovs = []
		self.maxcovs = []
		if n <= self.count:
			self.names = names
			self.rmsds = rmsds
			self.lengths = lengths
			self.mincovs = mincovs
			self.maxcovs = maxcovs
		else:
			lasti = None
			for i in np.arange(0, n, n/self.count):
				if int(i) == lasti: continue
				lasti = int(i)
				self.names.append(names[int(i)])
				self.rmsds.append(rmsds[int(i)])
				self.lengths.append(lengths[int(i)])
				self.mincovs.append(mincovs[int(i)])
				self.maxcovs.append(maxcovs[int(i)])

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

def plot_3d_densities(kde, fig, ax, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100):
	X, Y = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([X.ravel(), Y.ravel()])

	Z = np.reshape(kde.evaluate(positions).T, X.shape)
	surf = ax.plot_surface(X, Y, Z, cmap='magma', linewidth=0, antialiased=False)

	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	#ax.set_zlim(-1.0, 1.0)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')
	ax.set_zlabel('p')
	fig.colorbar(surf, ax=ax)

def plot_3d_posteriors(pkde, nkde, fig, ax, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100):
	X, Y = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([X.ravel(), Y.ravel()])

	posteriors = pkde.evaluate(positions) * PRIOR / (pkde.evaluate(positions) * PRIOR + nkde.evaluate(positions) * (1-PRIOR))

	Z = np.reshape(posteriors.T, X.shape)

	aspect = (rmsdlim[1] - rmsdlim[0]) / (mincovlim[1] - mincovlim[0]) / 1.414
	surf = ax.plot_surface(X, Y, Z, cmap='magma', linewidth=0, antialiased=False)
	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	ax.set_zlim(0.0, 1.0)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')
	ax.set_zlabel('p')
	ax.set_title('Posteriors (n={}+{})'.format(len(pkde.mincovs), len(nkde.mincovs)))
	fig.colorbar(surf, ax=ax)

def plot_posteriors(pkde, nkde, fig, ax, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100, colorbar=True):
	X, Y = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([X.ravel(), Y.ravel()])

	posteriors = pkde.evaluate(positions) * PRIOR / (pkde.evaluate(positions) * PRIOR + nkde.evaluate(positions) * (1-PRIOR))

	Z = np.reshape(posteriors.T, X.shape)

	aspect = (rmsdlim[1] - rmsdlim[0]) / (mincovlim[1] - mincovlim[0]) / 1.414
	im = ax.imshow(np.rot90(Z), cmap='magma', extent=rmsdlim+mincovlim, aspect=aspect)
	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')

	if colorbar: fig.colorbar(im, ax=ax)

	ax.set_title('Posteriors (n={}+{})'.format(len(pkde.mincovs), len(nkde.mincovs)))
	
def plot_independent_posteriors(positive, negative, figure, ax, rmsdlim, mincovlim, resolution=100, rbw=(1., 1.), cbw=(1., 1.)):
	rmsds, mincovs = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, \
		mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([rmsds.ravel(), mincovs.ravel()])

	pos_rkde = scipy.stats.gaussian_kde(positive.rmsds)
	pos_ckde = scipy.stats.gaussian_kde(positive.mincovs)
	neg_rkde = scipy.stats.gaussian_kde(negative.rmsds)
	neg_ckde = scipy.stats.gaussian_kde(negative.mincovs)
	pos_rkde.set_bandwidth(bw_method=pos_rkde.factor*rbw[0])
	pos_ckde.set_bandwidth(bw_method=pos_ckde.factor*cbw[0])
	neg_rkde.set_bandwidth(bw_method=neg_rkde.factor*rbw[1])
	neg_ckde.set_bandwidth(bw_method=neg_ckde.factor*cbw[1])

	rvalues = pos_rkde.evaluate(positions[0]) * PRIOR / (pos_rkde.evaluate(positions[0]) * PRIOR + neg_rkde.evaluate(positions[0]) * (1 - PRIOR))
	cvalues = pos_ckde.evaluate(positions[1]) * PRIOR / (pos_ckde.evaluate(positions[1]) * PRIOR + neg_ckde.evaluate(positions[1]) * (1 - PRIOR))
	values = rvalues * cvalues
	Z = np.reshape(values.T, rmsds.shape)

	aspect = (rmsdlim[1] - rmsdlim[0]) / (mincovlim[1] - mincovlim[0]) / 1.414
	im = ax.imshow(np.rot90(Z), cmap='magma', extent=rmsdlim+mincovlim, aspect=aspect)

	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')
	ax.set_title('Posteriors (n={}+{})'.format(len(positive.mincovs), len(negative.mincovs)))
	figure.colorbar(im, ax=ax)


def plot_3d_independent_posteriors(positive, negative, fig, ax, rmsdlim, mincovlim, resolution=200, rbw=(1., 1.), cbw=(1., 1.)):
	rmsds, mincovs = np.mgrid[rmsdlim[0]:rmsdlim[1]:resolution*1j, \
		mincovlim[0]:mincovlim[1]:resolution*1j]
	positions = np.vstack([rmsds.ravel(), mincovs.ravel()])

	pos_rkde = scipy.stats.gaussian_kde(positive.rmsds)
	pos_ckde = scipy.stats.gaussian_kde(positive.mincovs)
	neg_rkde = scipy.stats.gaussian_kde(negative.rmsds)
	neg_ckde = scipy.stats.gaussian_kde(negative.mincovs)
	pos_rkde.set_bandwidth(bw_method=pos_rkde.factor*rbw[0])
	pos_ckde.set_bandwidth(bw_method=pos_ckde.factor*cbw[0])
	neg_rkde.set_bandwidth(bw_method=neg_rkde.factor*rbw[1])
	neg_ckde.set_bandwidth(bw_method=neg_ckde.factor*cbw[1])

	rvalues = pos_rkde.evaluate(positions[0]) * PRIOR / (pos_rkde.evaluate(positions[0]) * PRIOR + neg_rkde.evaluate(positions[0]) * (1 - PRIOR))
	cvalues = pos_ckde.evaluate(positions[1]) * PRIOR / (pos_ckde.evaluate(positions[1]) * PRIOR + neg_ckde.evaluate(positions[1]) * (1 - PRIOR))
	values = rvalues * cvalues
	Z = np.reshape(values.T, rmsds.shape)

	surf = ax.plot_surface(rmsds, mincovs, Z, cmap='magma', linewidth=0, antialiased=False)
	ax.set_xlim(rmsdlim)
	ax.set_ylim(mincovlim)
	ax.set_zlim(0.0, 1.0)
	ax.set_xlabel('RMSD')
	ax.set_ylabel('Coverage')
	ax.set_zlabel('p')
	ax.set_title('Posteriors (n={}+{})'.format(len(positive.mincovs), len(negative.mincovs)))
	fig.colorbar(surf, ax=ax)

def plot_univariate_densities(dataset, fig, ax1, ax2, lim1, lim2, resolution=200, rbw=1., cbw=1.):
	X1 = np.arange(lim1[0], lim1[1], (lim1[1] - lim1[0])/resolution)
	X2 = np.arange(lim2[0], lim2[1], (lim2[1] - lim2[0])/resolution)

	rkde = scipy.stats.gaussian_kde(dataset.rmsds)
	ckde = scipy.stats.gaussian_kde(dataset.mincovs)
	#Y1 = rkde.evaluate(X1)
	#Y2 = ckde.evaluate(X2)

	#ax1.plot(X1, Y1)
	#ax2.plot(X2, Y2)

	rkde.set_bandwidth(bw_method=rkde.factor*rbw)
	Y1 = rkde.evaluate(X1)
	ax1.plot(X1, Y1)
	ckde.set_bandwidth(bw_method=ckde.factor*cbw)
	Y2 = ckde.evaluate(X2)
	ax2.plot(X2, Y2)

	ax1.set_xlim(lim1)
	ax2.set_xlim(lim2)
	

def plot_univariate_posteriors(pvalues, nvalues, fig, ax, lim, resolution=200, xlabel='', ylabel='', bw=(1., 1.)):
	X = np.arange(lim[0], lim[1], (lim[1]-lim[0])/resolution)
	pkernel = scipy.stats.gaussian_kde(pvalues)
	nkernel = scipy.stats.gaussian_kde(nvalues)

	pkernel.set_bandwidth(bw_method=pkernel.factor*bw[0])
	nkernel.set_bandwidth(bw_method=nkernel.factor*bw[1])
	posteriors = pkernel.evaluate(X)  * PRIOR / (pkernel.evaluate(X) * PRIOR + nkernel.evaluate(X) * (1-PRIOR))
	ax.plot(X, posteriors)
	ax.set_xlim(lim)
	ax.set_ylim([0, 1])
	ax.set_xlabel(xlabel)
	ax.set_ylabel(ylabel)
	#ax.plot(kde.rmsds, kde.mincovs, c=(1., 1., 1., 0.004), marker='.', linewidth=0)


#def plot_densities(positive, negative, unknown, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=100, density_outfile='kde_plot.png', scatter_outfile='kde_scatter.png', text='', posterior_outfile='kde_posterior.png', univar_post_outfile='univar_posterior.png', post_surf_outfile='posterior_surface.png', unklabel='Unknown', dens_surf_outfile='dens_surf_plot.png', univar_dens_outfile='univar_dens_plot.png', rbw=(1.75,1.75,1.75), cbw=(1.25,2.,1.5), indep_post_outfile='indep_post_plot.png', dpi=1200):
def plot_densities(positive, negative, unknown, rmsdlim=(0., 6.), mincovlim=(0., 1.), resolution=400, density_outfile='kde_plot.png', scatter_outfile='kde_scatter.png', text='', posterior_outfile='kde_posterior.png', univar_post_outfile='univar_posterior.png', post_surf_outfile='posterior_surface.png', unklabel='Unknown', dens_surf_outfile='dens_surf_plot.png', univar_dens_outfile='univar_dens_plot.png', rbw=(1.75,1.75,1.75), cbw=(1.25,2.,1.5), indep_post_outfile='indep_post_plot.png', dpi=600, univar_surf_outfile='univar_surf_plot.png', univar_postsurf_outfile='univar_postsurf.png'):

	print(rbw, cbw)
	##posteriors in 2D
	#figure = Figure()
	#canvas = FigureCanvas(figure)
	#ax = figure.add_subplot(1, 1, 1)
	#plot_posteriors(positive, negative, fig=figure, ax=ax, rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	#figure.savefig(posterior_outfile, dpi=dpi)
	#exit()

	#a bunch of posteriors in 2D
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []

	poscbws = (1., 2., 3., 4., 5., 6., 7., 8.)
	poscbws = np.arange(1., 7., 1.)
	negcbws = (1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13., 14.)
	negcbws = np.arange(1., 9., 1.)

	for i in range(len(poscbws)*len(negcbws)):
		ax.append(figure.add_subplot(len(poscbws), len(negcbws), i+1))
		ax[-1].tick_params(labelsize=3.)
	#figure.tight_layout()
	i = 0
	for poscbw in poscbws:
		positive.kernel.set_bandwidth(bw_method=positive.kernel.factor * poscbw)
		for negcbw in negcbws:
			negative.kernel.set_bandwidth(bw_method=negative.kernel.factor * negcbw)
			plot_posteriors(positive, negative, fig=figure, ax=ax[i], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution/2, colorbar=False)
			ax[i].set_title('{}, {}'.format(poscbw, negcbw), fontsize=4.)
			ax[i].set_xlabel('', fontsize=3.)
			ax[i].set_ylabel('', fontsize=3.)
			ax[i].plot([2.8], [0.68], marker='+', color=(0.,1.,0.,0.3))
			i += 1
			negative.kernel.set_bandwidth(bw_method=negative.kernel.factor / negcbw)
		positive.kernel.set_bandwidth(bw_method=positive.kernel.factor / poscbw)
	figure.savefig(posterior_outfile, dpi=dpi)

	#univariate densities
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []
	plots = []
	for i in range(6):
		ax.append(figure.add_subplot(2, 4, i+1))
		if i % 2:
			ax[-1].set_xlabel('Coverage', fontsize=6)
		else:
			ax[-1].set_xlabel('RMSD', fontsize=6)
		ax[-1].set_ylabel('p', fontsize=6)
		ax[-1].tick_params(labelsize=6)
	ax.append(figure.add_subplot(2, 2, 4))
	ax[-1].tick_params(labelsize=0)
	ax[-1].axis([0, 1, 0, 1])
	ax[-1].text(0.0, 0.0, text, fontsize=5., fontname='monospace')


	plot_univariate_densities(positive, figure, ax[0], ax[1], lim1=rmsdlim, lim2=mincovlim, resolution=2*resolution, rbw=rbw[0], cbw=cbw[0])
	ax[0].set_title('Positive (n={})'.format(len(positive.mincovs)))
	plot_univariate_densities(negative, figure, ax[2], ax[3], lim1=rmsdlim, lim2=mincovlim, resolution=2*resolution, rbw=rbw[1], cbw=cbw[1])
	ax[2].set_title('Negative (n={})'.format(len(negative.mincovs)))
	plot_univariate_densities(unknown, figure, ax[4], ax[5], lim1=rmsdlim, lim2=mincovlim, resolution=2*resolution, rbw=rbw[2], cbw=cbw[2])
	ax[4].set_title('{} (n={})'.format(unklabel, len(unknown.mincovs)))
	figure.savefig(univar_dens_outfile, dpi=dpi)

	##univariate posteriors
	#figure = Figure()
	#canvas = FigureCanvas(figure)
	#ax = []
	#plots = []
	#for i in range(2): 
	#	ax.append(figure.add_subplot(1, 2, len(ax)+1))
	#	plots.append(Plot(fig=figure, canvas=canvas, ax=ax[-1]))

	#plot_univariate_posteriors(positive.rmsds, negative.rmsds, figure, ax[0], lim=rmsdlim, resolution=2*resolution, xlabel='RMSD', ylabel='p', bw=rbw)
	#plot_univariate_posteriors(positive.mincovs, negative.mincovs, figure, ax[1], lim=mincovlim, resolution=2*resolution, xlabel='Coverage', ylabel='p', bw=cbw)
	#figure.savefig(univar_post_outfile, dpi=dpi)

	#a bunch of univariate posteriors
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []
	plots = []
	for i in range(2): 
		ax.append(figure.add_subplot(1, 2, len(ax)+1))
		plots.append(Plot(fig=figure, canvas=canvas, ax=ax[-1]))

	for pcbw in (1., 2., 3., 4., 5., 6.):
		for ncbw in (8., 9.):
			plot_univariate_posteriors(positive.rmsds, negative.rmsds, figure, ax[0], lim=rmsdlim, resolution=2*resolution, xlabel='RMSD', ylabel='p', bw=rbw)
			plot_univariate_posteriors(positive.mincovs, negative.mincovs, figure, ax[1], lim=mincovlim, resolution=2*resolution, xlabel='Coverage', ylabel='p', bw=(pcbw, ncbw))
			figure.savefig(univar_post_outfile, dpi=dpi)

	##2D univariate posterior plot
	#figure = Figure()
	#canvas = FigureCanvas(figure)
	#ax = figure.add_subplot(1, 1, 1)
	#plot_independent_posteriors(positive, negative, figure, ax, rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution, rbw=rbw, cbw=cbw)
	#figure.savefig(univar_postsurf_outfile, dpi=dpi)

	#density plots
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []
	plots = []
	for i in range(4): 
		ax.append(figure.add_subplot(2, 2, len(ax)+1))
		plots.append(Plot(fig=figure, canvas=canvas, ax=ax[-1]))
		#print(dir(ax[-1]))

	plot_kde(positive, fig=figure, ax=ax[0], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[0].set_title('Positive (n={})'.format(len(positive.mincovs)))
	plot_kde(negative, fig=figure, ax=ax[1], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[1].set_title('Negative (n={})'.format(len(negative.mincovs)))
	plot_kde(unknown, fig=figure, ax=ax[2], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution)
	ax[2].set_title('{} (n={})'.format(unklabel, len(unknown.mincovs)))

	ax[3].axis([0, 1, 0, 1])
	ax[3].text(0.0, 0.0, text, fontsize=5., fontname='monospace')

	figure.savefig(density_outfile, dpi=dpi)
	plot_rmsd_cov(positive, fig=figure, ax=ax[0])
	plot_rmsd_cov(negative, fig=figure, ax=ax[1])
	plot_rmsd_cov(unknown, fig=figure, ax=ax[2])
	figure.savefig(scatter_outfile, dpi=dpi)

	#posteriors in 3D
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = figure.add_subplot(1, 1, 1, projection='3d')
	plot_3d_posteriors(positive, negative, fig=figure, ax=ax, rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution*BEAUTYFACTOR)
	figure.savefig(post_surf_outfile, dpi=dpi)

	#densities in 3D
	figure = Figure()
	canvas = FigureCanvas(figure)
	ax = []
	for i in range(3):
		ax.append(figure.add_subplot(2, 2, i+1, projection='3d'))
	ax.append(figure.add_subplot(2, 2, 4))

	plot_3d_densities(positive, fig=figure, ax=ax[0], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution*BEAUTYFACTOR)
	ax[0].set_title('Positive (n={})'.format(len(positive.mincovs)))
	plot_3d_densities(negative, fig=figure, ax=ax[1], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution*BEAUTYFACTOR)
	ax[1].set_title('Negative (n={})'.format(len(negative.mincovs)))
	plot_3d_densities(unknown, fig=figure, ax=ax[2], rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=resolution*BEAUTYFACTOR)
	ax[2].set_title('{} (n={})'.format(unklabel, len(unknown.mincovs)))

	maxz = max(ax[0].get_zlim()[1], ax[1].get_zlim()[1], ax[2].get_zlim()[1])
	ax[0].set_zlim([ax[0].get_zlim()[0], maxz])
	ax[1].set_zlim([ax[1].get_zlim()[0], maxz])
	ax[2].set_zlim([ax[2].get_zlim()[0], maxz])
	ax[3].axis([0, 1, 0, 1])
	ax[3].text(0.0, 0.0, text, fontsize=5., fontname='monospace')
	figure.savefig(dens_surf_outfile, dpi=dpi)

	#3D univariate posterior plots
	#figure = Figure()
	#canvas = FigureCanvas(figure)
	#ax = figure.add_subplot(1, 1, 1, projection='3d')
	#plot_3d_independent_posteriors(positive, negative, figure, ax, rmsdlim=rmsdlim, mincovlim=mincovlim, resolution=BEAUTYFACTOR*resolution, rbw=rbw, cbw=cbw)
	#figure.savefig(univar_surf_outfile, dpi=dpi)


def plot_rmsd_cov(kde, fig, ax):
	ax.plot(kde.rmsds, kde.mincovs, c=(1., 1., 1., 0.01), marker='.', linewidth=0)


def main(positivefn, negativefn, unknownfn, count=1000, density_outfile='density_plot.png', scatter_outfile='scatter_plot.png', stretch=1, posterior_outfile='posterior_plot.png', univar_post_outfile='univar_posterior_plot.png', post_surf_outfile='post_surf_plot.png', min_present=50, unklabel='Unknown', dens_surf_outfile='dens_surf_plot.png', univar_dens_outfile='univar_dens_plot.png', indep_post_outfile='indep_post_plot.png', univar_surf_outfile='univar_surf_plot.png', univar_postsurf_outfile='univar_postsurf_plot.png'):
	with open(positivefn) as f: positive = Dataset(f, count=count, mode='any', marg=stretch, min_present=min_present)
	info('Positives: n = {}'.format(len(positive.mincovs)))
	positive.gen_rmsd_mincov_kde()

	for name, rmsd, mincov in zip(positive.names, positive.rmsds, positive.mincovs):
		#if mincov < 0.5: 
		if mincov < 0.5:
			print('{}\t{:0.2f}\t{:0.2%}'.format(name, rmsd, mincov))
	#exit()
	#with open(negativefn) as f: negative = Dataset(f, count=count, mode='any', marg=stretch, min_present=min_present)
	with open(negativefn) as f: negative = Dataset(f, count=count, mode='any', marg=stretch, min_present=min_present)
	info('Negatives: n = {}'.format(len(negative.mincovs)))
	negative.gen_rmsd_mincov_kde()
	with open(unknownfn) as f: unknown = Dataset(f, count=count, mode='any', marg=stretch, min_present=min_present)
	unknown.gen_rmsd_mincov_kde()
	info('Unknown: n = {}'.format(len(unknown.mincovs)))

	text = '#{}, {}, {}. stretch={}, n={}\n'.format(positivefn, negativefn, unknownfn, stretch, count)

	min_rmsd, max_rmsd = 0.0, 8.0
	min_mincov, max_mincov = 0.4, 1.0

	unkpoints = np.vstack([unknown.rmsds, unknown.mincovs])

	#this eliminates the low-cov peak
	#rbw = (1.75, 1.75, 1.75)
	#cbw = (1., 2., 1.5)

	#previous settings
	#cbw = (1.0, 2.25, 1.5)

	#this eliminates the low-cov peak, altering as little as possible
	
	rbw = (1.25, 1.25, 1.25)
	cbw = (2.0, 3.0, 1.5)

	#
	#rbw = (2.25, 2.25, 2.25)
	#cbw = (1., 2., 1.5)

	cbw = (1., 2., 1.5)

	text += '#Bandwidth multipliers: Pos {}, Neg {}\n'.format(cbw[0], cbw[1])
	positive.kernel.set_bandwidth(bw_method=positive.kernel.factor*cbw[0])
	negative.kernel.set_bandwidth(bw_method=positive.kernel.factor*cbw[1])
	unknown.kernel.set_bandwidth(bw_method=positive.kernel.factor*cbw[2])

	unkposteriors = positive.evaluate(unkpoints) * PRIOR / (positive.evaluate(unkpoints) * PRIOR + negative.evaluate(unkpoints) * (1 - PRIOR))
	#text += '#Alignment\tPosterior\tRMSD\tCoverage\n'
	text += '#Alignment'.ljust(32)
	text += 'Posterior'.ljust(12)
	text += 'RMSD'.ljust(6)
	text += 'Coverage\n'
	#for name, post, rmsd in zip(unknown.names, unkposteriors, unknown.rmsds): 
		#if post < 0.7: continue
		#print('{}\t{}\t{}'.format(name, post, rmsd))

	positive.kernel.set_bandwidth(bw_method=positive.kernel.factor/cbw[0])
	negative.kernel.set_bandwidth(bw_method=positive.kernel.factor/cbw[1])
	unknown.kernel.set_bandwidth(bw_method=positive.kernel.factor/cbw[2])

	for post, name, rmsd, mincov in sorted(zip(unkposteriors, unknown.names, unknown.rmsds, unknown.mincovs))[::-1][:10]:
		#text += '{}\t{:0.02e}\t{:0.1f}\t{:0.0%}\n'.format(name, post, rmsd, mincov)
		text += '{}'.format(name).ljust(32)
		text += '{:0.02e}'.format(post).ljust(12)
		text += '{:0.4f}'.format(rmsd).ljust(9)
		text += '{:0.0%}\n'.format(mincov)
	print(text)
	plot_densities(positive, negative, unknown, rmsdlim=(min_rmsd, max_rmsd), mincovlim=(0.0, max_mincov), resolution=100, 
		unklabel=unklabel, 
		density_outfile=density_outfile, 
		scatter_outfile=scatter_outfile, 
		posterior_outfile=posterior_outfile, 
			text=text, 
		univar_dens_outfile=univar_dens_outfile, 
		univar_post_outfile=univar_post_outfile, 
		post_surf_outfile=post_surf_outfile,
		dens_surf_outfile=dens_surf_outfile, 
		indep_post_outfile=indep_post_outfile,
		univar_surf_outfile=univar_surf_outfile,
		rbw=rbw, 
		cbw=cbw)


	#print('stats: ({} <= RMSD <= {}, {} <= minCov <= {}, n={})'.format(min_rmsd, max_rmsd, min_mincov, max_mincov, count*3))
	#print('positive:', positive.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))
	#print('negative:', negative.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))
	#print('unknown:', unknown.kernel.integrate_box([min_rmsd, min_mincov], [max_rmsd, max_mincov]))

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument('positive')
	parser.add_argument('negative')
	parser.add_argument('unknown')
	parser.add_argument('--min-present', type=int, default=50)
	parser.add_argument('-n', type=int, default=1000)
	parser.add_argument('-o', nargs=9, default=['density_plot.png', 'scatter_plot.png', 'posterior_plot.png', 'univar_post_plot.png', 'post_surf_plot.png', 'dens_surf_plot.png', 'univar_dens_plot.png', 'indep_post_plot.png', 'univar_surf_plot.png', 'univar_postsurf_plot.png'])
	parser.add_argument('-s', default=1, type=int)
	parser.add_argument('--unklabel', default='Unknown')
	parser.add_argument('--beauty-factor', type=int, default=BEAUTYFACTOR, help='extra resolution for 3d plots (warning: drastically increases runtime)')

	parser.add_argument('--prior', type=float)

	args = parser.parse_args()

	if args.prior is None: PRIOR = 0.5
	else: PRIOR = args.prior

	BEAUTYFACTOR = 8 if args.beauty_factor is None else args.beauty_factor

	main(args.positive, args.negative, args.unknown, count=args.n, density_outfile=args.o[0], stretch=args.s, scatter_outfile=args.o[1], posterior_outfile=args.o[2], univar_post_outfile=args.o[3], post_surf_outfile=args.o[4], dens_surf_outfile=args.o[5], min_present=args.min_present, unklabel=args.unklabel, univar_dens_outfile=args.o[6], indep_post_outfile=args.o[7], univar_surf_outfile=args.o[8], univar_postsurf_outfile=args.o[9])
