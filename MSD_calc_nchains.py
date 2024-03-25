#!/usr/bin/python

# Script:  msdions.py
# Purpose: calculate mean sq displacement of various types
# Syntax:  msdions.py < filename
# Example: msdions.py < test.dump (dump file with scaled coordinates)
# Author:  Lisa Hall from Mark's g(r) code

# derived from fortran code
# -------------------------------------------------------

import sys,string
from numpy import *
from math import *
import fileinput

# INPUT PARAMETERs
nconf =  20 #20030 # number of configurations to take data for
npoly =    2000  #need this for com calc

def msdcalc(): #defining the function to calculate msd for all types and chains
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6, msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly
    global natoms, typea,msdcall, chain, msdvec8, chain_id
    global x0_chain, y0_chain, z0_chain, x_chain, y_chain, z_chain
    # needs xc,yc,zc,xbox,ybox,zbox,drg
    # outputs g[]
    #from numpy import *
    # debug
##    r2mol=zeros(numatoms,float32) #todo:change to num molecules
##    ninmol=zeros(numatoms, int)
    msd[1]=0
    #print("what type",len(msd))
    msd[2]=0
    #print("what type 2",len(msd))
    msd[3]=0
    msd[4]=0
    msd[5]=0
    msd[6]=0
    msd[7]=0
    msd[8]=0
    """For chain"""
    for ii in range(1,chain):
        xr = x_chain[ii] - x0_chain[ii]
        yr = y_chain[ii] - y0_chain[ii]
        zr = z_chain[ii] - z0_chain[ii]
        r2 = xr*xr + yr*yr + zr*zr
        #print("what is this r2", r2)
        msd[8] = msd[8]+r2
        #print("msd8", msd[8])


    #msd[1]=msd[1]/num1
    #msd[2]=msd[2]/num2
    #msd[3]=msd[3]/num3
    #msd[4]=msd[4]/num4
    #msd[5]=msd[5]/num5
    #msd[6]=msd[6]/num6
    msd[7]=msd[7]/npoly
    msd[8]=msd[8]/chain
    msdvec1[msdcall]=msdvec1[msdcall]+msd[1] #msdvector, we can do it for a chain
    #print("what is msdvec1", (msdvec1))
    msdvec2[msdcall]=msdvec2[msdcall]+msd[2]
    #print("what is msdvec2", (msdvec2))
    msdvec3[msdcall]=msdvec3[msdcall]+msd[3]
    #print("what is msdvec3", (msdvec3))
    msdvec4[msdcall]=msdvec4[msdcall]+msd[4]
    msdvec5[msdcall]=msdvec5[msdcall]+msd[5]
    msdvec6[msdcall]=msdvec6[msdcall]+msd[6]
    msdvec7[msdcall]=msdvec7[msdcall]+msd[7]
    msdvec8[msdcall]=msdvec8[msdcall]+msd[8]
    #print("what is msdvec8", (msdvec8))
    timestep[msdcall]=step

    #OUT = open(file, 'a')
   # OUT.write("%7i %8.4f %8.4f %8.4f\n" % (step,msd[1],msd[2],msd[3]))
   # for ii in range(1,natoms+1):
    #    OUT.write("%7i %7i %7i %7i %7i %8.4f %8.4f %8.4f\n" % (mol[ii],typea[ii],xi[ii],yi[ii],zi[ii],xc[ii],yc[ii],zc[ii]))
   # OUT.close()
# end def msdcalc

def wholecalculation(): #function to read through the file
    global msdvec1, msdvec2, msdvec3, msdvec4, msdvec5, msdvec6,msdvec7, timestep, msd, xx, yy, zz, x0, y0, z0, x0com, y0com, z0com, xcom, ycom, zcom, step, num1, num2, num3, num4, num5, num6, npoly
    global skip, chain, msdvec8, x0_chain, y0_chain, z0_chain
    global natoms, typea,msdcall, x_chain, y_chain, z_chain, chain_id

    infiles = ["production.dump"] #PSsorted.lammpstrj
    file = 'msdcm_R1skip%d.csv'% skip
    dummy=0
    msdcall=0
    msdvec1=zeros(nconf,float)
    #print("what is this",msdvec1)
    #print("\n")
    msdvec2=zeros(nconf,float)
    msdvec3=zeros(nconf,float)
    msdvec4=zeros(nconf,float)
    msdvec5=zeros(nconf,float)
    msdvec6=zeros(nconf,float)
    msdvec7=zeros(nconf,float)
    msdvec8=zeros(nconf,float)
    timestep=zeros(nconf,int)

    IN = fileinput.input(infiles)
    natoms = 0
    num_counts = {}
    #Not sure if this relevant
    for loopnum in range(0,1):
        IN.readline()
        line = IN.readline()      # time step
        print("skip timestep %s" % line)
        IN.readline()
        line = IN.readline()      # number of atoms
        fields =  line.split() #string.split(line)
        natoms = int(fields[0])
        dim=natoms+1
        IN.readline() #information on box size
        line = IN.readline()
        line = IN.readline()
        line = IN.readline()
        line = IN.readline()
        #print()
        """We can write a statement that check the longest chain"""
        for j in range(1,dim):
            traj = IN.readline()

            [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(traj)
            chain_id = int(molj)
            if chain_id in num_counts:
                num_counts[chain_id] +=1
            elif chain_id not in num_counts:
                num_counts[chain_id]=1
        print("num_counts", (num_counts))
        print("\n")
        print("\n")
        print("num_counts", len(num_counts))
        print("\n")
        print("\n")

        # val_list = list(num_counts.values()) #chain lengths
        # print("the longest chain is ", max(val_list))
        # print("the val list is", val_list)
        # print("the length list is", len(val_list))
        # new_list = []
        # for i in val_list:
        #     if i >=123:
        #     #if 100 <= i <= 120:
        #         new_list.append(val_list.index(i)+1)
        #         continue
        # new_list
        # print("\n")
        # print("\n")
        # print("the new list is", new_list)
        # print("\n")
        # print("\n")
        # print("the new length is", len(new_list))
        mol = []
        for key, value in num_counts.items():
            #if  value == 140:
            if 130 <= value <= 150:
            #if  10 < value < 215:
                mol.append(key)

        print("list of IDs of chains with 140 beads is",mol)
        print("\n")
        print("\n")
        print("length of new", len(mol))
        #print("mol_id for", val_list.index(22)+1)
            #print("first line is", traj)
        for j in range(1,dim):
            line = IN.readline()
            #print("first line is", line)
        #continue
    chain = max(num_counts.values())
    chain_id = max(num_counts)


    for loopnum in range(0,1): #read starting timestep

        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        IN.readline()
        line = IN.readline()      # time step
        #print("next line", line)
        fields = string.split(line)
        step = int(fields[0])
        #print("read step %i" % step)
        IN.readline()
        line = IN.readline()      # number of atoms
        fields = string.split(line)
        natoms = int(fields[0])  #reads number of atoms from first configuration; don't support changing num atoms
        dim=natoms+1
        xc=zeros(dim,float32)
        yc=zeros(dim,float32)
        zc=zeros(dim,float32)
        xx=zeros(dim,float32)
        yy=zeros(dim,float32)
        zz=zeros(dim,float32)
        x0=zeros(dim,float32)
        y0=zeros(dim,float32)
        z0=zeros(dim,float32)
        x0com=zeros(npoly+1,float32)
        y0com=zeros(npoly+1,float32)
        z0com=zeros(npoly+1,float32)
        nbead=zeros(npoly+1)
        xi=zeros(dim)
        yi=zeros(dim)
        zi=zeros(dim)
        typea=[0]*dim
        chain_id = [0]*chain

        mol=[0]*dim #chain id
        msd=zeros(9,float32)
        IN.readline()
        line = IN.readline()
        [xm,xp] = list(map(float,line.split()))
        line = IN.readline()
        [ym,yp] = list(map(float,line.split()))
        line = IN.readline()
        [zm,zp] = list(map(float,line.split()))
        line = IN.readline()
        xbox = xp - xm
        ybox = yp - ym
        zbox = zp - zm
        vol = xbox*ybox*zbox
        xbox2 = xbox/2.0
        xcom=zeros(npoly+1,float32) #initialize these before every config read
        ycom=zeros(npoly+1,float32)
        zcom=zeros(npoly+1,float32)
        xcm=ycm=zcm=0.0
        x0_chain = y0_chain = z0_chain = zeros(chain, float32)
        x_chain = y_chain = z_chain = zeros(chain, float32)
        xcm_chain=ycm_chain=zcm_chain= 0
        num1 = 0
        num2 = 0
        num3 = 0
        num4 = 0
        num5 = 0
        num6 = 0
        for j in range(1,dim):
            line = IN.readline()
            #todo: make a list called "dumpstyle" so this can be changed easily in one place?
            [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)
            k=int(ii)
            typea[k] = int(typej)
            if typea[k] == 1 and molj == chain_id:
                num1=num1+1
            elif typea[k] == 2:
                num2=num2+1
            elif typea[k] == 3:
                num3=num3+1
            elif typea[k] == 4:
                num4=num4+1
            elif typea[k] == 5:
                num5=num5+1
            elif typea[k] == 6:
                num6=num6+1
            mol[k] = int(molj)
            xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
            yc[k] = ybox*(float(x2)-0.5)
            zc[k] = zbox*(float(x3)-0.5)
            xi[k] = int(round(float(n1)))
            yi[k] = int(round(float(n2)))
            zi[k] = int(round(float(n3)))
            xx[k] = xc[k] + xbox*xi[k]
            yy[k] = yc[k] + ybox*yi[k]
            #print("yy is",len(yy))
            zz[k] = zc[k] + zbox*zi[k]
        for m in range(1, chain):
            xcm_chain = xcm_chain + xx[m]
            #print("xcm_chain is ",  xcm_chain)
            ycm_chain = ycm_chain + yy[m]
            #print("ycm_chain is ",  ycm_chain)
            zcm_chain = zcm_chain + zz[m]
            #print("new length is", len(zz))
        xcm_chain = xcm_chain/chain
        #print("xcm_chain final", xcm_chain)
        ycm_chain = ycm_chain/chain
        #print("ycm_chain final", ycm_chain)
        zcm_chain = zcm_chain/chain

        for j in range(1, chain):
            x_chain[j] = xx[j] - xcm_chain
            y_chain[j] = yy[j] - ycm_chain
            z_chain[j] = zz[j] - zcm_chain
            #print("x_chain from center of mass is ",  x_chain)

        for i in range(1, chain):
            x0_chain[i] = xx[i]
            #print("ycm_chain final", x0_chain[i])
            y0_chain[i] = yy[i]
            z0_chain[i] = zz[i]




            istep=step

        msdcalc()
        msdcall+=1

        #print("Reading config file....")

        istart = loopnum + 1
        #print("loopnum", loopnum)
        # Read configuration from zconfig
        for kconf in range(istart,nconf):
          try:
            IN.readline()
            line = IN.readline()      # time step
            fields = string.split(line)
            step = int(fields[0])
            IN.readline()
            line = IN.readline()      # number of atoms
            fields = string.split(line)
            num = fields[0]
            IN.readline()
            line = IN.readline()
            [xm,xp] = list(map(float,line.split()))
            line = IN.readline()
            [ym,yp] = list(map(float,line.split()))
            line = IN.readline()
            [zm,zp] = list(map(float,line.split()))
            line = IN.readline()
            xbox = xp - xm
            ybox = yp - ym
            zbox = zp - zm
            vol = xbox*ybox*zbox
            xbox2 = xbox/2.0
            xcom=zeros(npoly+1,float32) #initialize these before every config read
            ycom=zeros(npoly+1,float32)
            zcom=zeros(npoly+1,float32)
            xcm=ycm=zcm=0.0
            nbead=zeros(npoly+1)
            for j in range(1,dim):
                line = IN.readline()
                #todo: make a list called "dumpstyle" so this can be changed easily in one place?
                [ii,molj,typej,q,x1,x2,x3,n1,n2,n3] = string.split(line)

                k=int(ii)
                typea[k] = int(typej)
                mol[k] = int(molj)
                if typea[k] == 1 and molj == chain_id:
                    xc[k] = xbox*(float(x1)-0.5) #scaled coords go from 0 to 1; want to go from -xbox/2 to xbox/2
                    yc[k] = ybox*(float(x2)-0.5)
                    zc[k] = zbox*(float(x3)-0.5)
                    xi[k] = int(round(float(n1)))
                    yi[k] = int(round(float(n2)))
                    zi[k] = int(round(float(n3)))
                    xx[k] = xc[k] + xbox*xi[k]
                    yy[k] = yc[k] + ybox*yi[k]
                    zz[k] = zc[k] + zbox*zi[k]

            for m in range(1, chain):
                xcm_chain = xcm_chain + xx[m]
                #print("xcm_chain is ",  xcm_chain)
                ycm_chain = ycm_chain + yy[m]
                zcm_chain = zcm_chain + zz[m]
                #print("new length is", len(zz))
            xcm_chain = xcm_chain/chain
            #print("xcm_chain final", xcm_chain)
            ycm_chain = ycm_chain/chain
            #print("ycm_chain final", ycm_chain)
            zcm_chain = zcm_chain/chain

            for j in range(1, chain):
                x_chain[j] = xx[j] - xcm_chain
                y_chain[j] = yy[j] - ycm_chain
                z_chain[j] = zz[j] - zcm_chain
            #print("x_chain from center of mass is ",  len(x_chain))



            # calculate msd for this config
            msdcalc()
            msdcall+=1
            #    Output
            #OUT = open(file, 'a')
            #OUT.write("%7i %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (step,msd[1],msd[2],msd[3],msd[4],msd[5]))
            #OUT.close()
            #except: break
            #print("%7i %7i %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f %8.4f\n" % (istep,step,msd[1],msd[2],msd[3],msd[4],msd[5],msd[6],msd[7]))
          except: break
        fileinput.close()

    # end of loop
    OUT = open(file, 'w')
    #OUT.write("msd\n")
    OUT.write("step, msd_type1, msd_type2, msd_type3, msd_type4, msd_type5, msd_type6, msd_com, msd_chain\n")
    for i in range(0,msdcall):
    ##    msdvec1[i]=msdvec1[i]/(nconf-i)
    ##    msdvec2[i]=msdvec2[i]/(nconf-i)
    ##    msdvec3[i]=msdvec3[i]/(nconf-i)
    ##    msdvec4[i]=msdvec4[i]/(nconf-i)
    ##    msdvec5[i]=msdvec5[i]/(nconf-i)
        OUT.write("%7i, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f, %8.4f\n" % (timestep[i],msdvec1[i],msdvec2[i],msdvec3[i],msdvec4[i],msdvec5[i],msdvec6[i],msdvec7[i], msdvec8[i]))
    OUT.close()

#main:
skip=0
wholecalculation()
