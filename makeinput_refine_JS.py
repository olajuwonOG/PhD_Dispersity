#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 10:19:18 2022

@author: taofeektejuosho
"""
#System 1 creation

import sys, string, math
import numpy as np
import random as rnd

#INPUT PARAMETERS
iseed = 500             # random number seed
nbeads = 1                # number of beads in monomer (3, 5, 7, or 9 for ionomer project)
#nxbeads = 3       #number of beads in the x segment of only backbone, to be randomly added in
#x_y = 0     #ratio of uncharged x segments per regular 'y' segments with charged monomers
nbeadsperpoly = 100
npoly = 1000
minsep = 1.0                # allowed separation in overlap check
#cisize = 0.5                #counterion diameter/bead diameter; to adjust density
#z_c = 1                   # counterion valence                 
dens = 0.85               # bead density

bond = 0.97  # bond length depends on bond potential, but close to 1 is good enough


# type names: charged monomer, neutral monomer, ion
#            0   1     2     3   4     5  6   7
atype = ( '00', 'NM', 'CM', 'IO', '00','00','00','XB' )

# define monomer

if (nbeads == 1): # This is type 1, meaning all monomers are the same = PE
    sequence = [1]
else:
    sys.exit ("define sequence of monomers.")

INPUT_LAMMPS = open('input.lammps', 'w')
INPUT_PDB = open('input_system_1.pdb', 'w')
INPUT_PSF = open('input_system_1.psf', 'w')

print(sequence)

nmonomers_ = []
for num1, num2 in zip(nbeadsperpoly, npoly):
    nmonomers_.append(num1*num2)
nmonomers = sum(nmonomers_)

print('Number of monomers:', nmonomers)

ntypes = 7 # Using type 1 in this case.

#Determining system size
ntot = (nmonomers)*nbeads #no counterions
vol = (ntot)/dens
side = vol**(1./3.)
print (side)
dim = ntot+1
# simulation cell parameters
hx = side
hy = side
hz = side

hx2 = hx/2.
hy2 = hy/2.
hz2 = hz/2.

vol = hx * hy * hz

nbonds = ntot-sum(npoly)

print('nbonds', nbonds)
print() 
print(("Total number of particles:",ntot))
print(("Number of chains =", npoly))
print(("beads in monomer =", nbeads))
print(("monomers total =", nmonomers))

print(("Number of atoms types = ",ntypes))
print(("seed = ", iseed))

print(" ")
print("Geometry:")
print(("dens = ", dens))

print(("vol = ", vol))

print(("metric: %10.4f %10.4f %10.4f\n\n" % (hx, hy, hz)))


# init position variables
xc=np.zeros((dim,),dtype=float) 
yc=np.zeros((dim,),dtype=float) 
zc=np.zeros((dim,),dtype=float) 
cx=np.zeros(dim)
cy=np.zeros(dim)
cz=np.zeros(dim)
print(xc)


# Build polymers #What does this mean?

rg2ave=0.0
rgave=0.0
rend2ave = 0.0
typeb=[0]*dim
molnum=[0]*dim
q=[0.0]*dim
i0=0
k=0
k1 = 0


for i in npoly:
    for ix in range(i):
        lengthcurrentpoly = 0
        for ii in nbeadsperpoly:
            for iy in range(ii):
                if npoly.index(i) == nbeadsperpoly.index(ii):
                    print(ii)
                    currentmonomer = ix*ii + iy
                    #print("current monomer is", currentmonomer)
                    seq=sequence
                    seqnum = 0
                    for iz in seq:
                        seqnum = seqnum + 1
                        k = k + 1
                        lengthcurrentpoly = lengthcurrentpoly + 1
                        typeb[k] = iz
                        print(typeb[k])
                        print(len(typeb))
                        molnum[k] = ix + 1
                        if iy == 0 and seqnum == 1:
                            k1 = k
                            xc[k] = rnd.random()*hx #num1 = random.randint(0, 9)
                            yc[k] = rnd.random()*hy
                            zc[k] = rnd.random()*hz
                        else:
                        # pick random direction; scale to be bond length 
                            dx = rnd.random()-0.5
                            dy = rnd.random()-0.5
                            dz = rnd.random()-0.5
                            r = np.sqrt(dx*dx+dy*dy+dz*dz)
                            scale = bond/r
                            dx = scale*dx
                            dy = scale*dy
                            dz = scale*dz

                            xc[k] = xc[k-1] + dx
                            yc[k] = yc[k-1] + dy
                            zc[k] = zc[k-1] + dz


  # calculate R and R_G
    # k2= k1 + lengthcurrentpoly-1
    # xcm = sum(xc[k1:k2+1])/lengthcurrentpoly
    # ycm = sum(yc[k1:k2+1])/lengthcurrentpoly
    # zcm = sum(zc[k1:k2+1])/lengthcurrentpoly
    # xg = xc[k1:k2+1]-xcm
    # yg = yc[k1:k2+1]-ycm
    # zg = zc[k1:k2+1]-zcm
    # rg2 = (np.dot(xg,xg) + np.dot(yg,yg) + np.dot(zg,zg))/lengthcurrentpoly
    #     # end to end
    # rend2 = (xc[k1]-xc[k2])**2 + (yc[k1]-yc[k2])**2 + (zc[k1]-zc[k2])**2
    # rend2ave = rend2ave + rend2
    # rg2ave = rg2ave + rg2
    # rgave = rgave + np.sqrt(rg2)
    # print ("current rg", rg2)
        

# for ix in range(npoly[1]):
#     lengthcurrentpoly = 0
#     for iy in range(nmonomersperpoly[1]):
#         currentmonomer = ix*nmonomersperpoly[1] + iy
#         # if currentmonomer in xmonomers:
#         #     seq=xsequence
#         # else:
#         seq=sequence
#         seqnum = 0
#         for iz in seq:
#             seqnum = seqnum + 1
#             k = k + 1
#             lengthcurrentpoly = lengthcurrentpoly + 1
#             typeb[k] = iz
#             molnum[k] = ix + 1
#             if iy == 0 and seqnum == 1:
#                k1 = k
#                xc[k] = rnd.random()*hx #num1 = random.randint(0, 9)
#                yc[k] = rnd.random()*hy
#                zc[k] = rnd.random()*hz
#             else:
#                 # pick random direction; scale to be bond length 
#                 dx = rnd.random()-0.5
#                 dy = rnd.random()-0.5
#                 dz = rnd.random()-0.5
#                 r = np.sqrt(dx*dx+dy*dy+dz*dz)
#                 scale = bond/r
#                 dx = scale*dx
#                 dy = scale*dy
#                 dz = scale*dz
                
#                 xc[k] = xc[k-1] + dx
#                 yc[k] = yc[k-1] + dy
#                 zc[k] = zc[k-1] + dz


#   # calculate R and R_G
#     k2= k1 + lengthcurrentpoly-1
#     xcm = sum(xc[k1:k2+1])/lengthcurrentpoly
#     ycm = sum(yc[k1:k2+1])/lengthcurrentpoly
#     zcm = sum(zc[k1:k2+1])/lengthcurrentpoly
#     xg = xc[k1:k2+1]-xcm
#     yg = yc[k1:k2+1]-ycm
#     zg = zc[k1:k2+1]-zcm
#     rg2 = (np.dot(xg,xg) + np.dot(yg,yg) + np.dot(zg,zg))/lengthcurrentpoly
#     # end to end
#     rend2 = (xc[k1]-xc[k2])**2 + (yc[k1]-yc[k2])**2 + (zc[k1]-zc[k2])**2
#     rend2ave = rend2ave + rend2
#     rg2ave = rg2ave + rg2
#     rgave = rgave + np.sqrt(rg2)
#     print ("current rg", rg2)    
    
    
    
    
# print("Polymers built.")
# rg2ave = rg2ave/npoly
# rgave = rgave/npoly
# rend2ave = rend2ave/npoly
# rave = rend2ave/rg2ave
# print(("<R_G^2> <R_G> = ",rg2ave,rgave))
# print(("<R_end^2>= ",rend2ave))


#no charges
# for k in range(1,ntot-ncounterions+1):
#     if (xc[k] > hx):
#         cx[k] = int(xc[k]/hx)
#         xc[k] = xc[k] - cx[k]*hx - hx2
#     elif (xc[k] < 0.0):
#         cx[k] = -int((-xc[k]+hx)/hx)
#         xc[k] = xc[k] - cx[k]*hx - hx2
#     else:
#         cx[k] = 0
#         xc[k] = xc[k] - hx2
#     if (yc[k] > hy):
#         cy[k] = int(yc[k]/hy)
#         yc[k] = yc[k] - cy[k]*hy - hy2
#     elif (yc[k] < 0.0):
#         cy[k] = -int((-yc[k]+hy)/hy)
#         yc[k] = yc[k] - cy[k]*hy - hy2
#     else:
#         cy[k] = 0
#         yc[k] = yc[k] - hy2
#     if (zc[k] > hz):
#         cz[k] = int(zc[k]/hz)
#         zc[k] = zc[k] - cz[k]*hz - hz2
#     elif (zc[k] < 0.0):
#         cz[k] = -int((-zc[k]+hz)/hz)
#         zc[k] = zc[k] - cz[k]*hz - hz2
#     else:
#         cz[k] = 0
#         zc[k] = zc[k] - hz2


INPUT_LAMMPS.write("#PolymerMelt Replica 3 - 11/2023\n")
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atoms\n" %     ntot)
INPUT_LAMMPS.write("%10i    bonds\n" %     nbonds)
INPUT_LAMMPS.write("%10i    angles\n" %     0)
INPUT_LAMMPS.write("%10i    dihedrals\n" % 0)
INPUT_LAMMPS.write("%10i    impropers\n" % 0)
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("%10i    atom types\n" % 1)
INPUT_LAMMPS.write("%10i    bond types\n" % 1)

INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write(" %16.8f %16.8f   xlo xhi\n" % (-hx2,hx2))
INPUT_LAMMPS.write(" %16.8f %16.8f   ylo yhi\n" % (-hy2,hy2))
INPUT_LAMMPS.write(" %16.8f %16.8f   zlo zhi\n" % (-hz2,hz2))
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Atoms\n")
INPUT_LAMMPS.write("\n")

mass = 1.0

# Atoms output
mass = 1.0

# Polymers
i = 0
imol = 0

for i in range(1,dim):
    itype = typeb[i]
    aname = atype[itype]

  # could use a dictionary here between type and segname
    if itype != 3:
        imol = molnum[i]
        segname = "POLY"
    elif itype == 3:
        #imol = npoly+1 #the molecule number for all counterions is the same; it's more like a group number
        imol = i-ntot+npoly #LMH now each ion has its own molecule number
        segname = "CION"

    INPUT_LAMMPS.write("%6i %6i %2i %6.2f %9.4f %9.4f %9.4f %6i %6i %6i\n" % (i, imol, typeb[i], q[i], xc[i], yc[i], zc[i], cx[i], cy[i], cz[i]))
   # INPUT_PDB.write("ATOM  %5i  %2s  NONE    1     %7.3f %7.3f %7.3f  1.00  0.00\n" %  (i,aname, xc[i], yc[i], zc[i] ))
   # INPUT_PSF.write("%8i %4s %3i  %2s   %2s   %2s   %8.6f       %7.4f %10i\n" %  (i,segname,imol,aname,aname,aname,q[i],typeb[i],0))

#Bonds   
INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Bonds\n")
INPUT_LAMMPS.write("\n")
jbond1 = np.zeros(nbonds+1)
jbond2 = np.zeros(nbonds+1)
pbond1 = np.zeros(nbonds+1)
pbond2 = np.zeros(nbonds+1)
ibond=0
pbond=0
i0 = 0
for i in range(1,ntot):
        #if not at the end of the polymer
        if molnum[i+1] == molnum[i]:
            ibond = ibond+1 #the bond number
            j=i+1
            INPUT_LAMMPS.write("%8i  1 %8i %8i\n" % (ibond,i,j))

INPUT_LAMMPS.write("\n")
INPUT_LAMMPS.write("Masses\n")
INPUT_LAMMPS.write("\n")

for ii in range(1,ntypes+1):
    INPUT_LAMMPS.write("%3i  1.0\n" % ii)

# #Close files
INPUT_LAMMPS.close()


print("LAMMPS, pdb, psf output complete.")
