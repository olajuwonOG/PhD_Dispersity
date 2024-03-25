#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:59:21 2023

@author: taofeektejuosho_che_desktop
"""
import resource
import sys
import numpy as np
import time
import fileinput
import pandas as pd
import csv
#Creating a function that calculates the average of the stress autocorrelation from the 3 components
stress = []
stress_new = []

with open("average.csv", 'r') as file:
#csvreader = csv.reader(file)
    for row in file:
        row_ = float(row.rstrip())
        stress_new.append((row_))

#stress_new = [float(i) for i in stress_]
#print(stress_new)


for i in range(1,len(stress_new)):
    low = int(round(i*.9))
    high = int(round(i*1.1))
    sum = 0
    tot = 0
    while(high > len(stress_new)):
        high -= 1
    for n in range(low,high):
        sum += stress_new[n-1]
        tot += 1
    if low != high:
        stress_new[i-1] = sum/tot
print(stress_new)
output = open("running_average3.csv", "w")
for i in stress_new:
#print(len(stress))
    output.write(str(i)+"\n")

print("code finished running")