#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 20:26:56 2023

@author: guterlj
"""

import numpy as np
from matplotlib.collections import PatchCollection

data = np.load("test_167196_bgp_mesh.npy",allow_pickle=True).tolist()
fig = plt.figure()
ax = fig.add_subplot(111)
patches = [matplotlib.patches.Polygon(np.array(v)) for v in data['mesh']]

p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4, edgecolor=None)
p.set_array(data['KTEBS'])
ax.add_collection(p)
ax.set_xlim([0,2])
ax.set_ylim([-3,3])

fig.colorbar(p)

fig.show()