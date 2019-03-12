#!/usr/bin/env python

from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

import sys, os

rmsdlim = [2.0, 3.5]
covlim = [0.8, 1.0]
qualim = [0.6, 0.8]

rmsdlim = [0.0, 6.0]
covlim = [0.0, 1.0]
qualim = [0.0, 1.0]

fig = plt.figure()
axs = []
axs.append(fig.add_subplot(221, projection='3d'))
axs.append(fig.add_subplot(222, projection='3d'))
axs.append(fig.add_subplot(223, projection='3d'))

ax4 = fig.add_subplot(224, projection='3d')

color = ['b', 'r', 'g']
color = [(0., 0., 1., 0.04), (1., 0., 0., 0.5), (0., 0.5, 0., 0.5), (1., 1., 0., 1.)]
i = 0
for ax, fn in zip(axs, sys.argv[1:]):

	labels = ['X', 'Y', 'Z']
	values = []
	with open(fn) as f:
		for l in f:
			if not l.strip(): 
				continue
			elif l.startswith('#'): 
				labels = l[1:].split('\t')
			else: 
				#print(l.split('\t'))
				row = [float(x) for x in l.split('\t')]
				if row[0] < 0: continue
				if row[0] < rmsdlim[0] or row[0] > rmsdlim[1]: continue
				if row[1] < covlim[0] or row[1] > covlim[1]: continue
				if row[2] < qualim[0] or row[2] > qualim[1]: continue
				values.append(row)

	values = np.array(values)

	ax.plot(values[:,0], values[:,1], values[:,2], linewidth=0, marker='.', color=color[i])
	ax4.plot(values[:,0], values[:,1], values[:,2], linewidth=0, marker='.', color=color[i])
	i += 1

	ax.set_xlabel(labels[0])
	ax.set_ylabel(labels[1])
	ax.set_zlabel(labels[2])

	ax.set_xlim(rmsdlim) #rmsd
	ax.set_ylim(covlim) #cov
	ax.set_zlim(qualim) #qual


	ax.set_title('{} n={:d}'.format(os.path.splitext(os.path.basename(fn))[0], len(values[:,0])))
ax4.set_xlabel(labels[0])
ax4.set_ylabel(labels[1])
ax4.set_zlabel(labels[2])

ax4.set_xlim(rmsdlim) #rmsd
ax4.set_ylim(covlim) #cov
ax4.set_zlim(qualim) #qual
if len(sys.argv) > 4: title = sys.argv[4]
else: title = ''
ax4.set_title(title)
plt.show()
