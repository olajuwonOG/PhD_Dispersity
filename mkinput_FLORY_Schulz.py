# Author: Taofeek Tejuosho

# lammps types
# 1 end cap
# 2 normal bead
# 3 side chain start
# 4 side chain
# 5 side chain end cap

import math
import random
import pandas as pd
from collections import Counter
import numpy as np
import scipy.stats as st

from scipy.stats import rv_discrete

def generateInputFile(nbeadsInput, monomers_per_polymer_input, npoly_input, minsep_input, dens_input, bond_length_input,
                      p_branch_end_input, filename, nmin, nmax, schulz_zimm=False, Flory_Schulz=True ):

    # class Flory_Schulz(rv_discrete):
    #       "Flory Schulz distribution"
    #       def _pmf(self, k, mu):
    #           return exp(-mu) * mu **k/factorial(k)
    # mu = 10
    # nmin = 10
    # nmax = 200


    class Flory_Schulz(st.rv_continuous):
         def _pdf(self, x, p):

    #         "Flory Schulz Distribution"
             return ((p**(x-1))*x*(1-p)**2)
             #return (x)*(p**2)*((1-p)**(x-1))
    p = 0.831
    iseed = 15
    my_cv = Flory_Schulz(a=nmin, b=nmax, name='Flory_Schulz')


    N = 1000
    monomers_per_polymer = my_cv.rvs(size = N, random_state = 1, p = p)
    polydispersity = True
    monomers_per_polymer = monomers_per_polymer.astype(int)
    list_mon = list(monomers_per_polymer)
    print(list_mon)
    print("length of the polymer", len(list_mon))
    print("Total number of beads",(monomers_per_polymer))
    """Calculating the dispersity"""
    z = Counter(list_mon)
    data = pd.DataFrame({"No_chains":z})
    print("This is the data ",data)
    data.index
    data.reset_index(inplace = True)
    data = data.rename(columns = {'index': "No_of_monomers"})
    data['weight_fraction'] = data["No_of_monomers"]*data["No_chains"]/sum(data["No_of_monomers"]* data["No_chains"])
    data['number_fraction'] = data["No_chains"]/sum(data["No_chains"])
    print(data)
    #display(data.to_string())
    with open("composition.csv", 'w') as file:
        #file.write("Index, No_of_monomers,  No_chains,  weight_fraction,  number_fraction")
        #file.write("\n")
        for i in data.index:
            NumberOfMonomers = str(data['No_of_monomers'][i])
            NumberFraction = str(data["number_fraction"][i])
            WeightFraction = str(data["weight_fraction"][i])
            file.write(NumberOfMonomers + "," + NumberFraction + "," + WeightFraction)
            file.write("\n")
    Mw = sum((data['No_of_monomers']*data['weight_fraction']))
    print("Average molecular weight is ", Mw)
    Mn = sum((data['No_of_monomers']*data['No_chains']))/sum(data['No_chains'])
    print("number average molecular weight is ", Mn)
    D = Mw/Mn
    print("Dispersity is ", D)

    pass
    # elif shape != 0:
        # iseed = 20
        # monomers_per_polymer = np.random.gamma(shape, scale=scale, size=npoly_input)
        # polydispersity = True
        # monomers_per_polymer = monomers_per_polymer.astype(int)
        # print("mon per pol", monomers_per_polymer.size)
        #for i in range(monomers_per_polymer.size):
        #    monomers_per_polymer[i] = int(monomers_per_polymer[i])

    # need more input verification

    if nbeadsInput > 0:
        print("input params detected")
        nbeads = nbeadsInput
        if not polydispersity:
            monomers_per_polymer = monomers_per_polymer_input  # monomers in polymer (including x)
        npoly = npoly_input  # number of polymers
        minsep = minsep_input  # minimum seperation to prevent overlap
        # pendant_size = 1.0  # pendant group diameter/bead diameter
        dens = dens_input  # bead density
        bond_length = bond_length_input  # bond length
        p_branch_end = p_branch_end_input
    else:
        print("no input params or incorrect params (7 fields)")
        print(
            "fields are:\nbeads per monomer (mainchain only)\nmonomers per polymer\nnumber of polymers\nminimum "
            "seperation to prevent overlap\nbead density\nbond length\nprobability of end of side branch")
        nbeads = 3  # beads per monomer mainchain
        monomers_per_polymer = 3  # monomers in polymer (including x)
        #npoly = 3  # number of polymers
        minsep = 1.0  # minimum seperation to prevent overlap
        # pendant_size = 1.0  # pendant group diameter/bead diameter
        dens = 0.85  # bead density
        bond_length = 0.97  # bond length
        p_branch_end = 0.5
    minsep2 = minsep ** 2
    #print(minsep2)
    #print("beads per monomer: ", nbeads)
    #print("number of beads in x segement of only backbone: ", nxbeads)
    #print("monomers per polymer: ", monomers_per_polymer)
    #print("Maximum chain length is", max(monomers_per_polymer))
    #print("Minimum chain length is", min(monomers_per_polymer))
    #print("number of polymers: ", npoly)
    #print("minimum seperation ", minsep)
    # print("pendant size: ", pendant_size)
    #print("bead density: ", dens)
    #random.seed(random_seed)

    #                        0                   1                   2                 3                4
    type_names = ['Polymer Endcap Bead', 'Main Chain Bead', 'Branch Start Bead', 'Branch Bead', 'Branch Endcap']
    bond_type_names = ['Endcap to Main Bead', 'Main bead to main bead',
                       'main bead to side chain start', 'side chain start to side chain middle',
                       'side chain middle to side chain end', 'side chain middle to side chain end']
    # branch cant occur at first or last main chain beads in the nbead monomer
    # for branching negative binomial distribution? Trials needed to get endcap defines final monomer
    bead_count = 0
    side_chain_total_length = 0
    polymers = []
    for i in range(npoly):
        monomers = []
        for j in range(monomers_per_polymer[i]):
            monomer = []
            # starting bead here (0)
            if j == 0:
                #TODO change starting bead back to endcap
                monomer.append(0)
                bead_count = bead_count + 1
            elif j == monomers_per_polymer[i] - 1:
                #TODO change back to endcap
                monomer.append(0)
                bead_count = bead_count + 1
            # regular main chain bead if not starting bead
            else:
                monomer.append(1)
                bead_count = bead_count + 1

            #print("what is this monomer", monomer)

            # append individual monomer to monomers list
            monomers.append(monomer)
        # append monomers that constitute polymer into polymer list
        polymers.append(monomers)
        #print("this is the final polymer", len(polymers))
    INPUT_LAMMPS = open(filename, 'w')
    if polydispersity:
        nmonomers = np.sum(monomers_per_polymer)
    else:
        nmonomers = monomers_per_polymer * npoly
    ntypes = len(type_names)
    # is this right?
    nbranches = nbeads * nmonomers
    ntot = bead_count
    vol = ntot / dens
    side = math.pow(vol, 1 / 3)
    dim = ntot
    hx, hy, hz = (side, side, side)
    hx2, hy2, hz2 = (hx / 2, hy / 2, hz / 2)
    if not polydispersity:
        nbonds = side_chain_total_length + (((monomers_per_polymer * nbeads) - 1) * npoly)
    else:
        ## TODO get working for branches
        nbonds = np.sum(monomers_per_polymer) - monomers_per_polymer.size
    """print()
    print("Total number of particles:", ntot)
    print("Number of chains =", npoly)
    print("actual density", ntot / vol)
    print("Average beads in main chain =", np.average(monomers_per_polymer) * nbeads)
    print("beads in side chain =", side_chain_total_length)
    print("average side chain length =", side_chain_total_length / nmonomers)
    print("monomers total =", nmonomers)
    print("number of bonds =", nbonds)
    print("Number of atoms types = ", ntypes)
    print("seed = ", random_seed)

    print()
    print("Geometry:")
    print("dens = ", dens)

    print("vol = ", vol)

    print()
    print("metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz))"""
    print("monomers total =", nmonomers)
    print("number of bonds =", nbonds)

    # position variables
    xc = np.zeros(dim)
    yc = np.zeros(dim)
    zc = np.zeros(dim)
    # velocity variables
    cx = np.zeros(dim)
    cy = np.zeros(dim)
    cz = np.zeros(dim)
    # Building polymers

    Rg2Ave = 0.0
    RgAve = 0.0
    Rend2Ave = 0.0
    types = np.zeros(dim)
    molnum = np.zeros(dim)
    q = np.zeros(dim)
    k = 0
    bonds = []

    for polymerNum in range(len(polymers)):
        lengthcurrentpoly = 0
        firstInPolymer = True
        for monomerNum in range(len(polymers[polymerNum])):
            inSideChain = False
            sideChainStart = (0, 0, 0)
            for unitNum in range(len(polymers[polymerNum][monomerNum])):
                unit = polymers[polymerNum][monomerNum][unitNum]
                types[k] = unit
                molnum[k] = polymerNum + 1
                if firstInPolymer:
                    xc[k] = random.random() * hx
                    yc[k] = random.random() * hy
                    zc[k] = random.random() * hz
                    firstInPolymer = False
                else:
                    if unitNum > 0:
                        prevUnit = polymers[polymerNum][monomerNum][unitNum - 1]
                    else:
                        prevUnit = polymers[polymerNum][monomerNum - 1][-1]
                    dx = random.random() - 0.5
                    dy = random.random() - 0.5
                    dz = random.random() - 0.5
                    r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
                    scale = bond_length / r
                    dx = scale * dx
                    dy = scale * dy
                    dz = scale * dz
                    if types[k] == 2:
                        inSideChain = True
                        xc[k] = xc[k - 1] + dx
                        yc[k] = yc[k - 1] + dy
                        zc[k] = zc[k - 1] + dz
                        sideChainStart = (xc[k], yc[k], zc[k])
                    elif types[k] == 4 or types[k] == 3:
                        xc[k] = xc[k - 1] + dx
                        yc[k] = yc[k - 1] + dy
                        zc[k] = zc[k - 1] + dz
                    elif (types[k] == 1 or types[k] == 0) and inSideChain:
                        xc[k] = sideChainStart[0] + dx
                        yc[k] = sideChainStart[1] + dy
                        zc[k] = sideChainStart[2] + dz
                        inSideChain = False
                    else:
                        xc[k] = xc[k - 1] + dx
                        yc[k] = yc[k - 1] + dy
                        zc[k] = zc[k - 1] + dz
                k = k + 1

    for k in range(types.size):
        if xc[k] > hx:
            cx[k] = int(xc[k] / hx)
            xc[k] = xc[k] - cx[k] * hx - hx2
        elif (xc[k] < 0.0):
            cx[k] = -int((-xc[k] + hx) / hx)
            xc[k] = xc[k] - cx[k] * hx - hx2
        else:
            cx[k] = 0
            xc[k] = xc[k] - hx2
        if yc[k] > hy:
            cy[k] = int(yc[k] / hy)
            yc[k] = yc[k] - cy[k] * hy - hy2
        elif yc[k] < 0.0:
            cy[k] = -int((-yc[k] + hy) / hy)
            yc[k] = yc[k] - cy[k] * hy - hy2
        else:
            cy[k] = 0
            yc[k] = yc[k] - hy2
        if zc[k] > hz:
            cz[k] = int(zc[k] / hz)
            zc[k] = zc[k] - cz[k] * hz - hz2
        elif zc[k] < 0.0:
            cz[k] = -int((-zc[k] + hz) / hz)
            zc[k] = zc[k] - cz[k] * hz - hz2
        else:
            cz[k] = 0
            zc[k] = zc[k] - hz2

    print("Polymers built.")

    INPUT_LAMMPS.write("#Polymer October/2022\n")
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("%10i    atoms\n" % ntot)
    INPUT_LAMMPS.write("%10i    bonds\n" % nbonds)
    # INPUT_LAMMPS.write("%10i    angles\n" %     0)
    # INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
    # INPUT_LAMMPS.write("%10i    impropers\n" % npendant)
    # INPUT_LAMMPS.write("%10i    impropers\n" % 0)
    INPUT_LAMMPS.write("\n")
    #TODO change back to full types
    INPUT_LAMMPS.write("1    atom types\n")
    INPUT_LAMMPS.write("1    bond types\n")
    #INPUT_LAMMPS.write("%10i    atom types\n" % ntypes)
    #INPUT_LAMMPS.write("%10i    bond types\n" % 2)
    ##INPUT_LAMMPS.write("%10i    improper types\n" % 1)

    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2, hx2))
    INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2, hy2))
    INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2, hz2))
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Atoms\n")
    INPUT_LAMMPS.write("\n")

    mass = 1.0

    for i in range(dim):
        unit = types[i]
        # could use a dictionary here between type and segname
        imol = molnum[i]  # this implies the polymers must be placed in before the counterions
        #TODO changed the default type to 1 (third input)

        INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (
            i+1, imol, 1, q[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))
    ##    INPUT_PDB.write("ATOM  %5i  %2s  NONE    1     %7.3f %7.3f %7.3f  1.00  0.00\n" %  (i,aname, xc[i], yc[i], zc[i] ))
    ##    INPUT_PSF.write("%8i %4s %3i  %2s   %2s   %2s   %8.6f       %7.4f %10i\n" %  (i,segname,imol,aname,aname,aname,q[i],typeb[i],0))

    # Bonds
    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Bonds\n")
    INPUT_LAMMPS.write("\n")
    bonds = np.zeros(nbonds)

    inSideChain = False
    bond_count = 0
    sideChainStart = 0
    for i in range(cx.size - 1):
        if molnum[i + 1] == molnum[i]:
            bond_count = bond_count + 1
            if (types[i + 1] == 1 and types[i] == 0) or (types[i + 1] == 0 and types[i] == 1):
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif (types[i+1]==0 and types[i]==0):
                bond_count = bond_count-1
            elif (types[i + 1] == 2 and types[i] == 0):
                inSideChain = True
                sideChainStart = i + 1
                INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif (types[i + 1] == 0) and inSideChain:
                inSideChain = False
                INPUT_LAMMPS.write("%8i  2 %8i %8i\n" % (bond_count, sideChainStart+1, i + 2))
            elif (types[i + 1] == 1 and types[i] == 1):
                #TODO change back to 3
                INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif types[i + 1] == 2 and types[i] == 1:
                inSideChain = True
                sideChainStart = i + 1
                INPUT_LAMMPS.write("%8i  4 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif (types[i + 1] == 1) and inSideChain:
                inSideChain = False
                INPUT_LAMMPS.write("%8i  4 %8i %8i\n" % (bond_count, sideChainStart+1, i + 2))
            elif types[i + 1] == 3 and types[i] == 2:
                INPUT_LAMMPS.write("%8i  5 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif types[i + 1] == 4 and types[i] == 2:
                INPUT_LAMMPS.write("%8i  6 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif (types[i + 1] == 3 and types[i] == 3):
                INPUT_LAMMPS.write("%8i  7 %8i %8i\n" % (bond_count, i+1, i + 2))
            elif types[i + 1] == 4 and types[i] == 3:
                INPUT_LAMMPS.write("%8i  8 %8i %8i\n" % (bond_count, i+1, i + 2))

    INPUT_LAMMPS.write("\n")
    INPUT_LAMMPS.write("Masses\n")
    INPUT_LAMMPS.write("\n")
    #TODO changed to default 1 type, 1 mass
    #for ii in range(ntypes):
    #    INPUT_LAMMPS.write("%3i  1.0\n" % ii)
    INPUT_LAMMPS.write("1  1.0\n")
    INPUT_LAMMPS.close()
    print("LAMMPS output complete.")


generateInputFile(1,0,100,1,0.85,0.97,0,'input.lammps',nmin = 3, nmax = 179)
