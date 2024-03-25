#!/usr/bin/env python
# coding: utf-8

# In[1]:


import math
import random
import pandas as pd
import numpy as np
import scipy.stats as st
from scipy.special import gamma, factorial
from scipy.stats import rv_continuous
from collections import Counter
from scipy.integrate import quad
import matplotlib.pyplot as plt


# # k.Mn

# In[33]:


def SZ(x,z,Mn):
        return ((z**(z+1)/gamma(z+1))*(x**(z)/Mn**(z+1))*np.exp(-z*x/Mn))

class Schulz_zimm(st.rv_continuous):
    def _pdf(self,x):
        self.z = 0.2
        self.Mn = 10
        self.integral, _ = quad(SZ, self.a, self.b, args=(self.z, self.Mn))
        z = 0.2
        Mn = 10
        return ((z**(z+1)/gamma(z+1))*(x**(z)/Mn**(z+1))*np.exp(-z*x/Mn))/self.integral
       
        #sigma = Mn/z , n = z
np.random.seed(10)
nmin = 3
nmax = 500
N = 1000
my_cv = Schulz_zimm(a=nmin, b=nmax, name='Schulz_Zimm')
print(my_cv)

monomers = my_cv.rvs(size = N)

#plt.hist(monomers, bins = 30)
#plt.axvline(monomers.mean(), color='k', linestyle='dashed', linewidth=2)
#plt.show()

monomers_per_polymer = monomers.astype(int)
list_mon = list(monomers_per_polymer)
print("length of the polymer", len(list_mon))
print(sum(monomers_per_polymer))
"""Calculating the dispersity"""
z = Counter(list_mon)
data = pd.DataFrame({"No_chains":z})
data.index
data.reset_index(inplace = True)
data = data.rename(columns = {'index': "No_of_monomers"})
data['weight_fraction'] = data["No_of_monomers"]*data["No_chains"]/sum(data["No_of_monomers"]* data["No_chains"])  
Mw = sum((data['No_of_monomers']*data['weight_fraction']))
print("Average molecular weight is ", Mw)
Mn = sum((data['No_of_monomers']*data['No_chains']))/sum(data['No_chains'])
print("number average molecular weight is ", Mn)
D = Mw/Mn
print("Dispersity is ", D)
a = np.array(list_mon)
# Creating histogram
fig, ax = plt.subplots(figsize =(10, 7))
ax.hist(a, bins = 60, edgecolor = "black", density = True)
#plt.title('SZ (MW = 80, D = 41.4)', fontsize = 30, fontweight = 'bold', fontname="Arial")
plt.xlabel("Number of Beads", fontsize = 30, fontname="Arial")
plt.ylabel("Probability", fontsize = 30, fontname="Arial")
plt.xticks(fontsize=20)
plt.yticks(fontsize=20)
plt.show()
# In[36]:


#samples = my_cv.rvs(size = 500)
#display(samples)


# In[37]:




