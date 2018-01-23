#!/usr/bin/env python2

"""
Hydration extraction (or any other cloud extraction) largely inspired by ANGULA from Luis Carlos Pardo:
https://gcm.upc.edu/en/members/luis-carlos/angula/ANGULA

MDAnalysis crucial function are largely presented here:
https://pythonhosted.org/MDAnalysis/documentation_pages/analysis/align.html

This is a short tested example.

Author: Mateusz Bieniek, bieniekmat@gmail.com
"""

import MDAnalysis
import numpy as np
from MDAnalysis.analysis.align import rotation_matrix


def coord_axis(origin, x_direction, third_atom):
    def normalise(vector):
        return vector / np.linalg.norm(vector)

    x_dim = normalise(origin - x_direction)
    tmp_vec = normalise(origin - third_atom)
    # cross the vectors
    y_dim = normalise(np.cross(x_dim, tmp_vec))
    z_dim = normalise(np.cross(x_dim, y_dim))

    # coordination axis (x, y, z) vectors normalised
    return x_dim, y_dim, z_dim


def pbc_correct(coord, pbc):
    """
    Move the coordinate to the right PBC.
    Assumes that the translation to the origin has been done.
    And that the selection is smaller than half the PBC.
    :param coord: coordinate along the pbc axis
    :param pbc: pbc size
    """
    if coord > pbc / 2:
        return coord - pbc
    elif coord < -(pbc / 2):
        return coord + pbc

    return coord

# input files
gro_path = 'my.gro'
xtc_path = 'my.xtc'
# output file
out_filepath = 'within5_OW_testpbc_originshifting.pdb'

# create the 'universe'
u = MDAnalysis.Universe(gro_path, xtc_path)

# open a 'writer' as W for the output
with MDAnalysis.Writer(out_filepath) as W:
    # extract all the IDs of HSS molecules
    hss_ids = set(u.select_atoms('resname HSS').atoms.resids)

    for hss_id in hss_ids:
        # for each time frame
        for ts in u.trajectory:
            print 'Time (ps)', ts.time

            sel_dst = 5   # in angstroms

            # e.g. around 5 (resid 1 and name S) -> 5 angstroms away from the Sulphur atom of the resid 1
            ow_and_hss = u.select_atoms('(name OW and around %d (resid %d and name S)) or resid %d'
                                        % (sel_dst, hss_id, hss_id))
            s_pos = ow_and_hss.atoms.S.position

            # move all the atoms to have Sulphur at (0, 0, 0)
            ow_and_hss.atoms.translate(-s_pos)

            # PBC correction
            for atom in ow_and_hss:
                corrected_x = pbc_correct(atom.position[0], u.dimensions[0])
                corrected_y = pbc_correct(atom.position[1], u.dimensions[1])
                corrected_z = pbc_correct(atom.position[2], u.dimensions[2])

                atom.position = (corrected_x, corrected_y, corrected_z)

            # HS1 is the X axis, S is the origin,
            sulphur_s_hs1_hs2_positions = ow_and_hss.select_atoms('name S', 'name HS1', 'name HS2').positions
            coord_axis_originS = coord_axis(*sulphur_s_hs1_hs2_positions)

            # reference frame: x, y z
            frame_ref = ((1, 0, 0), (0, 1, 0), (0, 0, 1))

            R, rmsd = rotation_matrix(coord_axis_originS, frame_ref)
            ow_and_hss.atoms.rotate(R)

            # output each .pdb frame
            W.write(ow_and_hss)
