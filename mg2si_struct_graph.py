#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 18 15:39:31 2022

@author: rlg3


Atomic Simulation Environment

Conda environment: alignn
"""

import ase.io as io
from ase.build import cut
from ase.spacegroup import crystal
from ase.data import colors
from jarvis.core.atoms import Atoms
from jarvis.core.graphs import Graph
import networkx as nx
import matplotlib.pyplot as plt
import dgl

a = 6.36
xtal = crystal(('Mg', 'Si'),
                       basis=[(0.25, 0.25, 0.25), (0, 0, 0)],
                       spacegroup=225,
                       cellpar=[a, a, a, 90, 90, 90])

#Change colors
colors.jmol_colors[12] = (189. / 255, 157. / 255, 227. / 255) # Mg
colors.jmol_colors[14] = (114. / 255, 173. / 255, 87. / 255)
colors.cpk_colors[12] = (189. / 255, 157. / 255, 227. / 255) # Mg
colors.cpk_colors[14] = (114. / 255, 173. / 255, 87. / 255)

# Create a new atoms instance with Co at origo including all atoms on the
# surface of the unit cell
mg2si = cut(xtal, origo=(0., 0., 0.), extend=1.01)

# Define the atomic bonds to show
bondatoms = []
symbols = mg2si.get_chemical_symbols()
for i in range(len(mg2si)):
    for j in range(i):
        if (mg2si.get_distance(i, j) < 3.2 and symbols[i] != symbols[j]):
            bondatoms.append((i, j))
        # elif (symbols[i] == symbols[j] == 'Si' and
        #       mg2si.get_distance(i, j) < 3.2):
        #     bondatoms.append((i, j))

# Create nice-looking image using povray
renderer = io.write('spacegroup-mg2si.pov', mg2si,
                    rotation='75y, 10x',
                    radii=0.4,
                    povray_settings=dict(transparent=False,
                                         camera_type='perspective',
                                         canvas_width=320,
                                         bondlinewidth=0.07,
                                         bondatoms=bondatoms))
#Mg2Si : FCC
atoms_dict = {'lattice_mat': [[3.891767547488381,
   8.710926e-10,
   2.2469128302502366],
  [1.2972555171714013, 3.66919375143634, 2.246912830250237],
  [9.018215e-10, 6.376841e-10, 4.493826657376472]],
 'coords': [[3.8917725, 2.7518925, 6.740737500000001],
  [1.2972575000000002, 0.9172975, 2.2469125],
  [0.0, 0.0, 0.0]],
 'elements': ['Mg', 'Mg', 'Si'],
 'abc': [4.493827, 4.493823, 4.49383],
 'angles': [60.0, 60.0001, 60.0],
 'cartesian': True,
 'props': ['', '', '']}

atoms = Atoms.from_dict(atoms_dict)


graph_fcc = Graph.atom_dgl_multigraph(
            atoms,
            cutoff=100.0,
            atom_features="atomic_number",
            max_neighbors=4,
            compute_line_graph=True,
            use_canonize=False,
        )


#Ti : HCP
atoms = Atoms.from_poscar('output_files/POSCAR_Ti')

graph_hcp = Graph.atom_dgl_multigraph(
            atoms,
            cutoff=100.0,
            atom_features="atomic_number",
            max_neighbors=12,
            compute_line_graph=True,
            use_canonize=False,
        )

G0 = dgl.to_networkx(graph_fcc[0])
nx.draw(G0, node_color = [(189. / 255, 157. / 255, 227. / 255),\
                          (189. / 255, 157. / 255, 227. / 255),\
                            (114. / 255, 173. / 255, 87. / 255)], arrows = False)

plt.savefig('graph0_mg2si.pdf', bbox_inches = 'tight')

plt.figure()
G1 = dgl.to_networkx(graph_fcc[1])
nx.draw(G1, node_color = 'xkcd:silver', edge_color = 'xkcd:silver', alpha = 0.3, arrows = False)

plt.savefig('graph1_mg2si.pdf', bbox_inches = 'tight', transparent = True)