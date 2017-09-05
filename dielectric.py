#!/usr/bin/env python
from astropy.io import ascii
import numpy as np
import sys #For input from bash shell
#print sys.argv[1]
import json
#Get the effective radius of the protein
filename = sys.argv[1]
file = ascii.read(filename)
i=0
length = len(file)
radgyr = [file[i][1] for i in range(0,length)]
effrad = np.power(1.29099*np.mean(radgyr),3)

#This block is to calculate the constant value used in Pitera 2001 or Simonson 1996
#x represents the denominator and y is for multiplying by 1/x

x = ((8.854e-12)*(1e-10)*(4*3.1415)*(1.381e-23)*(300))/((1.6e-19)**2)
y = 1/x


#Calculate the dipole moment of the protein
filename = sys.argv[2]
coords = ascii.read(filename)
dipole = [[coords[i][1+j]-coords[i][4+j] for j in range(0,3)] for i in range(0,len(coords))]
accumulated_variance = [0, 0 ,0]
half_var=[0,0,0]
half_dc=0
dielectric_constant = 0
total_variance = 0
half_total_var=0
x_axis=[coords[i][1]-coords[i][4] for i in range(0,len(coords))]
y_axis=[coords[i][2]-coords[i][5] for i in range(0,len(coords))]
z_axis=[coords[i][3]-coords[i][6] for i in range(0,len(coords))]

half_x=[coords[i][1]-coords[i][4] for i in range(0,len(coords)/2)]
half_y=[coords[i][2]-coords[i][5] for i in range(0,len(coords)/2)]
half_z=[coords[i][3]-coords[i][6] for i in range(0,len(coords)/2)]


#X coordinate
half_var[0]=np.var(half_x)
accumulated_variance[0] = np.var(x_axis)

#Y coordinate
half_var[1]=np.var(half_y)
accumulated_variance[1] = np.var(y_axis)

#Z coordinate
half_var[2]=np.var(half_z)
accumulated_variance[2] = np.var(z_axis)

half_total_var=half_var[0]+half_var[1]+half_var[2]
half_dc=( 161 + (( 89920*half_total_var)/effrad))/( 161 - (( 562*half_total_var )/effrad))
total_variance = accumulated_variance[0] + accumulated_variance[1] + accumulated_variance[2]
dielectric_constant = ( 161 + (( 89920*total_variance)/effrad))/( 161 - (( 562*total_variance )/effrad))



print half_dc
print dielectric_constant

f = open('dc.txt','w')
f.write('dielectric constant')
f.write('	')
json.dump( format(dielectric_constant, '.1f'), f )
f.close()










