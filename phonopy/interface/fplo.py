# Copyright (C) 2016 Marc Höppner
# All rights reserved.
#
# This file is part of phonopy.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in
#   the documentation and/or other materials provided with the
#   distribution.
#
# * Neither the name of the phonopy project nor the names of its
#   contributors may be used to endorse or promote products derived
#   from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
# FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
# COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import sys
import numpy as np

 
from phonopy.file_IO import iter_collect_forces 
from phonopy.interface.vasp import get_scaled_positions_lines, get_drift_forces
from phonopy.units import Bohr
from phonopy.cui.settings import fracval
from phonopy.structure.atoms import Atoms, symbol_map
from phonopy.structure.symmetry import Symmetry 
from phonopy.file_IO import parse_disp_yaml

def parse_set_of_forces(num_atoms,
                        forces_filenames,
                        verbose=False):
    hook = 'CONVERGED'
    force_sets = []
    for fplo_filename in forces_filenames:
        fplo_forces = iter_collect_forces(fplo_filename,
                                          num_atoms,
                                          hook,
                                          [6, 7, 8],
                                          word='TOTAL NEUTRAL FORCES')

        # rotate forces from the FPLO xyz coordinate system
        # into the coordinate system of phonopy

        # read the phonopy lattice vectors
        conf, cell = parse_disp_yaml(return_cell=True)
        ppy_lat = cell.get_cell()

        # Get the transformation matrix.
        #symmetry = Symmetry(cell)
        #T = symmetry.get_dataset()["transformation_matrix"]

        # read the fplo output file
        f = open(fplo_filename, 'r')
        outf = f.readlines()
        f.close
 
        # search lattice vectors
        fp_lat = np.zeros((3,3))
        for j in np.arange(len(outf)):
            if outf[j].find("lattice vectors") == 0:
                for i in np.arange(3):
                    fp_lat[i,:] = list(map(float, outf[j+i+1].split(":")[1].split())) # [bohr]

        if np.allclose(fp_lat, np.zeros((3,3))):
            sys.exit("(EE) The lattice vectors were not found in: " + fplo_filename)

        T = np.dot(np.linalg.inv(ppy_lat), fp_lat)

        if not np.allclose(np.dot(T.T, T), np.eye(3)):
            sys.exit("(EE) The basis transformation is wrong!" )

        for i in np.arange(len(fplo_forces)):
            fplo_forces[i] = np.dot(T, np.array(fplo_forces[i]))

        # Calculate the drift_force (sum over all forces has to be zero) and
        # normalize the forces.
        drift_force = get_drift_forces(fplo_forces, fplo_filename, verbose)
        force_sets.append(np.array(fplo_forces) - drift_force)

    return force_sets

def read_fplo(filename):
    from phonopy.structure.cells import get_cell_matrix 

    fplo_in = open(filename).readlines()

    nsort = None

    for i in np.arange(len(fplo_in)):
        # spacegroup
        if fplo_in[i].find("spacegroup")+1:
            curl_pos = fplo_in[i+1].find("{")+1
            spg = int(fplo_in[i+1][curl_pos:].split(",")[0])     # integer number of the spacegroup

            if spg != 1:
                sys.exit("(EE) Only SPG 1 / P1 implemented so far.")
        
        # unit of length
        if fplo_in[i].find("lengthunit")+1:
            curl_pos = fplo_in[i+1].find("{")+1
            if int(fplo_in[i+1][curl_pos:].split(",")[0]) == 1:
                lfactor = 1.0                                        # bohr
            elif int(fplo_in[i+1][curl_pos:].split(",")[0]) == 2:    
                lfactor = 1.0/Bohr                                   # Angstroem
            else:
                sys.exit("(EE) Length unit unknown: " + fplo_in[i+1])
         
        # lattice constants
        if fplo_in[i].find("lattice_constants")+1:
            curl_pos = fplo_in[i].find("{")+1
            a, b, c  = map(float, fplo_in[i][curl_pos:].split(","))
       
        # lattice angles
        if fplo_in[i].find("axis_angles")+1:
            curl_pos           = fplo_in[i].find("{")+1
            alpha, beta, gamma = map(float, fplo_in[i][curl_pos:].split(","))
          
        # number of Wyckoff positions       
        if fplo_in[i].find("nsort=")+1:
            nsort = int(fplo_in[i][:-2].split("=")[1])

        # Wyckoff positions
        if fplo_in[i].find("wyckoff_positions")+1:
            if nsort == None: sys.exit("(EE) nsort not found!")

            species   = np.zeros((nsort), dtype=int)
            positions = np.zeros((nsort, 3), dtype=float)

            for j in np.arange(nsort):
                ind  = i + j + 2
                indc = fplo_in[ind].find(",")
                indb = fplo_in[ind].find("}")

                species[j]     = symbol_map[fplo_in[ind][indc-3:indc-1]]
                positions[j,:] = list(map(float, fplo_in[ind][indc+2:indb].split(",")))
     
    lattice = get_cell_matrix(a, b, c, alpha, beta, gamma)
    cell = Atoms(numbers=species,
                 cell=lattice,
                 scaled_positions=positions)

#    symmetry = Symmetry(cell)
#    cell = Atoms(numbers=symmetry.get_dataset()["std_types"],
#                 cell=symmetry.get_dataset()["std_lattice"],
#                 scaled_positions=symmetry.get_dataset()["std_positions"])

    return cell, False


def write_fplo(directory, cell, sym=None):
    from phonopy.structure.cells import get_angles, get_cell_parameters 
     
    FEDIT='fedit14.00-47-x86_64'
    FPLO='fplo14.00-47-x86_64'

    import os

    lattice = cell.get_cell()
    positions = cell.get_scaled_positions()
    chemical_symbols = cell.get_chemical_symbols()

    alpha, beta, gamma = get_angles(lattice)
    a, b, c = get_cell_parameters(lattice)
 
    os.mkdir(directory)
    os.system('cp =.in ' + directory + '/')
    os.chdir(directory)

    # Write the structure plan to the pipe file.
    fpipe = open("=.t_pipe", "w")
    fpipe.write("@+@\n")

    # add spacegroup
    if sym != None:
        fpipe.write("@s@\n")
        fpipe.write("@{0:d}@\n".format(sym["spg"]))
        fpipe.write("@x@\n")

    # change lattice units to bohr
    fpipe.write("@u@\n")
    fpipe.write("@b@\n")
    fpipe.write("@x@\n")

    fpipe.write("@l@ {0:+10.6f} {1:+10.6f} {2:+10.6f}\n".format(a,b,c))
    fpipe.write("@a@ {0:+8.4f} {1:+8.4f} {2:+8.4f}\n".format(alpha, beta, gamma))
    fpipe.write("@n@ {0:d}\n".format(len(positions)))
    for i in np.arange(len(positions)):
        fpipe.write("@{0:d}@ {1} @ {2:16.8f} {3:16.8f} {4:16.8f}\n".format(i+1, chemical_symbols[i],
                                                                           positions[i,0], positions[i,1], positions[i,2]))
    fpipe.write("@+@\n")
    fpipe.write("@x@\n")
    fpipe.write("@q@\n")
    fpipe.write("@y@")
    fpipe.close()

    # Execute fedit with the pip file.
    os.system(FEDIT+ ' -p '+FPLO+' -pipe < ./=.t_pipe 2>./+log 1>/dev/null')

    # Clean up.
    os.remove("=.t_pipe")
    os.chdir("../")


def write_supercells_with_displacements(supercell,
                                        cells_with_displacements,
                                        sym=False):
    # write plain supercell
    write_fplo("supercell", supercell)

    if sym == False:
        for i, cell in enumerate(cells_with_displacements):
            write_fplo("supercell-%03d" % (i + 1), cell)

    else:
        print("Symmetrized supercells:")

        for i, cell in enumerate(cells_with_displacements):
            # find spacegroup
            sym      = dict()
            symmetry = Symmetry(cell, 1e-5)
            sym["spg"] = symmetry.get_dataset()["number"]
            print("  Supercell: {0:03d}".format(i+1))
            print("    Spacegroup: " + symmetry.get_international_table())
            print("    Number of non-equivalent atoms in %s-%03d/=.in: %d" % (
                "supercell", i + 1, len(symmetry.get_independent_atoms())))


            print(cell.get_cell())
            print(symmetry.get_dataset()["std_lattice"])

            # create cell with symmetry
            ineq = np.array(symmetry.get_independent_atoms())
            cell = Atoms(numbers=symmetry.get_dataset()["std_types"][ineq],
                         cell=symmetry.get_dataset()["std_lattice"],
                         scaled_positions=symmetry.get_dataset()["std_positions"][ineq])
            write_fplo("supercell-%03d" % (i + 1), cell, sym)
