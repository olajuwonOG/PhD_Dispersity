""" Radius of Gyration
    Garrett Levine
    June 11, 2014
    python radiusofgyration.py < equil.dump
    This program will compute the mean square and mean radius of gyration for each polymer."""

import sys,string
from numpy import *
from math import *

nconf = 200 #number of configurations to average over
iskip = 0
npoly = 1000
nbeads = 60 #this is the number of beads in a polymer chain


#Initial conditions
conf = 0
natoms = npoly*nbeads
avgrg = [0]*npoly
avgsquaredrg = [0]*npoly
totavgsquaredrg = 0
totavgrg = 0
stdev1 = 0
stdev2 = 0

file = 'RadGyration.txt'

def radgy(): #a function to calculate radius of gyration
    for jj in range(0,npoly):
        #print("first jj is", jj)
        for kk in range(1,nbeads+1): #calculating center of mass - Rcm
           k = nbeads*jj + kk
           zcm[jj] += zc[k]/nbeads
           #print("zc is the, ", zc[jj])
           ycm[jj] += yc[k]/nbeads
           xcm[jj] += xc[k]/nbeads
    #print(zcm[0])
    #raise Exception("stop")
    for jj in range(0,npoly): #Going through each chain
        for kk in range(1, nbeads+1): #Going through the beads in a chain to calculate R - Rcm
            k = nbeads*jj + kk
            rg[jj] += (((xc[k]-xcm[jj])**2)+((yc[k]-ycm[jj])**2)+((zc[k]-zcm[jj])**2))/nbeads
            squarerootrg[jj] = sqrt(rg[jj])
        #print("rg is ", rg[jj])


for i in range(0,iskip+1): #reading through the dump file
    sys.stdin.readline()
    sys.stdin.readline() #timestep
    sys.stdin.readline()
    line = sys.stdin.readline() # number of atoms
    fields =  line.split() #split the readline named line into element separated by spaces
    natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
    dim=natoms+1
    xc=zeros(dim,float32)
    yc=zeros(dim,float32)
    zc=zeros(dim,float32)
    zci=zeros(dim,float32)
    cx=zeros(dim)
    cy=zeros(dim)
    cz=zeros(dim)
    typea=[0]*dim
    mol=[0]*dim
    chain=[0]*dim
    imf=[0]*dim
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    rg=[0]*npoly
    squarerootrg = [0]*npoly
    sys.stdin.readline()
    line = sys.stdin.readline() #xbox bounds
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() #ybox bounds
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() #zbox bounds
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = str.split(line) #atom number, molecule number, bead type, charges, x,y,z, image flags

    conf = conf + 1
    print(j, "configuration is ", conf)

OUT = open(file, 'w')

print("reading config file...")
istart = iskip+1

for i in range(istart,nconf+iskip+1):
    sys.stdin.readline()
    sys.stdin.readline() #timestep
    sys.stdin.readline()
    line = sys.stdin.readline() # number of atoms
    fields = str.split(line) #reads the next line
    natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
    dim=natoms+1
    xc=zeros(dim,float32)
    yc=zeros(dim,float32)
    zc=zeros(dim,float32)
    zci=zeros(dim,float32)
    cx=zeros(dim)
    cy=zeros(dim)
    cz=zeros(dim)
    typea=[0]*dim
    mol=[0]*dim
    chain=[0]*dim
    imf=[0]*dim
    zcm=[0]*npoly
    xcm=[0]*npoly
    ycm=[0]*npoly
    rg=[0]*npoly
    squarerootrg = [0]*npoly
    settings = sys.stdin.readline()
    line = sys.stdin.readline() #xbox bounds
    [xm,xp] = map(float,line.split())
    line = sys.stdin.readline() #ybox bounds
    [ym,yp] = map(float,line.split())
    line = sys.stdin.readline() #zbox bounds
    [zm,zp] = map(float,line.split())
    line = sys.stdin.readline()
    xbox = xp - xm
    ybox = yp - ym
    zbox = zp - zm
    vol = xbox*ybox*zbox
    xbox2 = xbox/2.0
    ybox2 = ybox/2.0
    zbox2 = zbox/2.0
    for j in range(1,dim):
        line = sys.stdin.readline()
        [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = str.split(line) #atom number, molecule number, bead type, x,y,z, image flags
        k=int(ii) #sets the index k equal to the atom/bead number
        chain[k] = molj #sets the chain number for each bead equal to molecule number
        q1=chain[k] #sets the index q equal to the chain number
        chain[k] = int(chain[k])
        xc[k] = xbox*float(x1)
        yc[k] = ybox*float(x2)
        zc[k] = zbox*float(x3)
        """By default, atom coords are written in a scaled format (from 0 to 1). That is, an
        value of 0.25 means the atom is at a location 1/4 of the distance from xlo to xhi of the box boundaries."""
        if xc[k] > xbox:              #xc[k] in dump file is not from 0 to 1, sometimes.
            xc[k] = xc[k] - xbox
        elif xc[k] < 0:
            xc[k] = xc[k] + xbox
        if yc[k] > ybox:              #yc[k] in dump file is not from 0 to 1, sometimes.
            yc[k] = yc[k] - ybox
        elif yc[k] < 0:
            yc[k] = yc[k] + ybox
        if zc[k] > zbox:              #zc[k] in dump file is not from 0 to 1, sometimes.
            zc[k] = zc[k] - zbox
        elif zc[k] < 0:
            zc[k] = zc[k] + zbox
        xc[k] = xbox*float(x1) + int(n1)*xbox
        yc[k] = ybox*float(x2) + int(n2)*ybox
        zc[k] = zbox*float(x3) + int(n3)*zbox
    #print("xc first is", xc)
            #should have assigned each bead a chain number
    radgy()
    conf = conf + 1
    print("what is configuration ",conf)

    for jj in range(0,npoly): #Average of each chain relating to a configuration
        avgsquaredrg[jj] += rg[jj]/nconf
        avgrg[jj] += squarerootrg[jj]/nconf

OUT = open(file, 'w')
OUT.write("#%7i\n" % (conf))
OUT.write("polymer number    rg squared    mean rg\n")
for jj in range(0, npoly): #Average of the whole chain
    totavgsquaredrg += avgsquaredrg[jj]/npoly
    totavgrg += avgrg[jj]/npoly
    OUT.write("%8.4f   %8.4f   %8.4f\n" % (jj+1, avgsquaredrg[jj], avgrg[jj]))
#stdev1 = std(avgsquaredrg)*1.96/sqrt(npoly)
#stdev2 = std(avgrg)*1.96/sqrt(npoly)
OUT.write("Total Average")
OUT.write("%8.4f %8.4f %8.4f %8.4f\n" % (totavgsquaredrg, totavgrg, stdev1, stdev2))
OUT.close()
