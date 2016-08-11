#!/usr/bin/env python
"""
Remove the water molecules which have at least one atom in the in Z region specified.

e.g. The following will remove all SOL solvent molecules which have at least one atom in the Z dimension between coordinates 2 and 6
remove_solatoms.py -f file.gro -o output.gro -b 2 -e 6
"""

import argparse
import sys
import os
import MDAnalysis

if __name__ == "__main__":
    # argument parser
    argparser = argparse.ArgumentParser(description="Removes the water/sol if any of the atoms reside withing the request Z dimension")
    argparser.add_argument("-f", help="filename .gro from which to remove", metavar="source .gro filename", required=True, type=str)
    argparser.add_argument("-o", help="name of the filename", metavar="output filename", required=True)
    argparser.add_argument("-b", help="start/begin of Z dimension in angstroms", metavar="float", required=True, type=float)
    argparser.add_argument("-e", help="end of Z dimension in angstroms", metavar="float", type=float, required=True)

    # parse arguments
    args = argparser.parse_args(sys.argv[1:])
    source_filename = args.f
    assert os.path.exists(source_filename), "Argument -f must point to a file"
    output_filename = args.o
    assert '.' in output_filename, 'The output file name has to have an extension that MDAnalysis is capabable of writing, including .gro or .pdb'
    zbegin = args.b
    zend = args.e
    assert zbegin < zend, "The first Z dimension value has to be smaller than the second"

    print 'Remove water molecules if any of their atoms reside in Z dim ', zbegin, zend

    u = MDAnalysis.Universe(source_filename)

    # solvent molecules with at least one atom in the Z interval to be excluded
    exile_list_resids = set()
    for atom in u.atoms:
        # is solvent?
        if atom.resname != 'SOL':
            continue

        zpos = atom.position[2]
        if zpos > zbegin and zpos < zend:
            exile_list_resids.add(atom.resid)

    if len(exile_list_resids) == 0:
        print 'No molecules to be removed'
        exit(0)

    print 'Number of SOL molecules to be removed: ', len(exile_list_resids)

    # atom indices without the excluded molecules
    survivors_ids = [atom.id for atom in u.atoms if atom.resid not in exile_list_resids]
    filtered = u.atoms[survivors_ids]

    if filtered.n_atoms == 0:
        print 'No atoms survived the SOL molecules removal. Exiting.'
        exit(0)

    filtered.write(output_filename)

