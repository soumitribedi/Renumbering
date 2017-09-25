#renumber a molecule according to the numbering of the reference molecule
#All files should be in cartesian coordinates and in xyz format
#first filename after scriptname in command line is the file that is to be renumbered
#second filename is the reference molecule

from scipy.optimize import fsolve
import numpy as np
import sys
import os
import re
import math

from sys import argv
file_ref = argv[1]      #first file from command line
file_2 = argv[2]        #second file from command line

def orientation(filename):             #function orientation will reorient the input file in a certain plane
	file = open(filename)
	print  filename
	lines = file.readlines()
	file.close()

	for j in range(len(lines)):
		natom = len(lines)-2   #natom gives the total number of atoms in xyz format file
		break
	print natom
	xcoordinate = []
	ycoordinate = []
	zcoordinate = []
	atom_type = []
	dist_origin = []
	new_atom_type = []
	new_xcoordinate = []
	new_ycoordinate = []
	new_zcoordinate = []
	rnew_atom_type = []
	rnew_xcoordinate = []
	rnew_ycoordinate = []
	rnew_zcoordinate = []
	r2new_atom_type = []
	r2new_xcoordinate = []
	r2new_ycoordinate = []
	r2new_zcoordinate = []
	r3new_atom_type = []
	r3new_xcoordinate = []
	r3new_ycoordinate = []
	r3new_zcoordinate = []
	atom_num = []

	dist = [[0 for i in range(natom)] for j in range(natom)]
	connect = [[0 for i in range(natom)] for j in range(natom)]
	
	for i,line in enumerate(lines):
		if i>1 : 
			(n,x,y,z) = line.split()
			X = float(x)
			Y = float(y)
			Z = float(z)
			xcoordinate.append(X)        #storing all the x coordinate values from the file in one list
			ycoordinate.append(Y)        #storing all the y coordinate values from the file in one list
			zcoordinate.append(Z)        #storing all the z coordinate values from the file in one list
			atom_type.append(n)          #storing all the atom types in one list

	for i in range(0,natom) :
		for j in range(0,natom) :
			dist_x = xcoordinate[i] - xcoordinate[j]
			dist_y = ycoordinate[i] - ycoordinate[j]
			dist_z = zcoordinate[i] - zcoordinate[j]
			dist[i][j]   = math.sqrt(dist_x**2+dist_y**2+dist_z**2)               #distance matrix
			#dist[i][j]   = "%.5f" % math.sqrt(dist_x**2+dist_y**2+dist_z**2)     #distance matrix contains distances upto 5 points of decimal
	p=0
	p2=[]
	for a in range(0,natom) :
		for b in range(0,natom) :
			if (0.8 < dist[a][b] < 1.6) &  (b != a):
				connect[a][b] = atom_type[a]+atom_type[b]           #connectivity matrix e.g., it will show 'NC' when nitrogen is connected to carbon

	for a in range(0,natom):
		for b in range(0,natom):
			if connect[a][b]== 'NC' :
				p=p+1
		p2.insert(a,p)                    #p2 is a matrix which shows the number of carbon atoms connected to the nitrogen atoms
		p=0

	for a in range(0,natom):
		if p2[a] == 3 :
			print a+1,atom_type[a], p2[a]
			first_x = xcoordinate[a]          #coordinates of the nitrogen that is connected to 3 carbon atoms. In this case it is a unique nitrogen
			first_y = ycoordinate[a]
			first_z = zcoordinate[a]
			atom1=a
	for a in range(0,natom) :
		new_xcoordinate.append(xcoordinate[a] - first_x)   #shifting the origin to the coordinates of this unique nitrogen
		new_ycoordinate.append(ycoordinate[a] - first_y)
		new_zcoordinate.append(zcoordinate[a] - first_z)
	def rotation_z(theta):            #function which returns the value of x coordinate after rotation by theta along z axis
		f1 = x1*math.cos(math.radians(theta)) - y1*math.sin(math.radians(theta))  
		return (f1)
	def rotation_x(beta):             #function which returns the value of z coordinate after rotation by theta along x axis
		f2 = y2*math.sin(math.radians(beta)) + z2*math.cos(math.radians(beta))  
		return (f2)
	def rotation_y(ceta):             #function which returns the value of x coordinate after rotation by theta along y axis
		f3 = z3*math.sin(math.radians(ceta)) + x3*math.cos(math.radians(ceta))  
		return (f3)

	for a in range(0,natom):
		if p2[a] == 2 :                 #a is the index of the unique nitrogen connected to 2 carbon atoms. We are rotating the molecule to bring this nitrogen in the yz plane
			print a+1,atom_type[a], p2[a]
			x1=new_xcoordinate[a]
			y1=new_ycoordinate[a]
			z1=new_zcoordinate[a]
			theta = fsolve(rotation_z,(0))   #built-in function to solve for the equation inside rotation_z for the root which is close to zero
			atom2=a

	y_atom2 = new_xcoordinate[atom2]*math.sin(math.radians(theta)) + new_ycoordinate[atom2]*math.cos(math.radians(theta))    #y coordinate after rotation along z axis
	if y_atom2 < 0:      #if rotation has taken place such that the nitrogen is in the negative yz plane
		theta = theta+180   #rotate in the opposite direction so that the nitrogen is in the positive plane

	for a in range(0,natom):
		rnew_xcoordinate.append(new_xcoordinate[a]*math.cos(math.radians(theta)) - new_ycoordinate[a]*math.sin(math.radians(theta)))    #new coordinates of the molecule after rotation
		rnew_ycoordinate.append(new_xcoordinate[a]*math.sin(math.radians(theta)) + new_ycoordinate[a]*math.cos(math.radians(theta)))
		rnew_zcoordinate.append(new_zcoordinate[a]*1 )
#now we are rotating the molecule such that this nitrogen sits on the y axis, so making the z coordinate of the atom zero
	x2=rnew_xcoordinate[atom2]
	y2=rnew_ycoordinate[atom2]
	z2=rnew_zcoordinate[atom2]
	beta = fsolve(rotation_x,(0))    

	for a in range(0,natom):
		r2new_xcoordinate.append(rnew_xcoordinate[a]*1 )          #new coordinates. Now the 2nd nitrogen atom sits on the y axis
		r2new_ycoordinate.append(rnew_ycoordinate[a]*math.cos(math.radians(beta)) - rnew_zcoordinate[a]*math.sin(math.radians(beta)))
		r2new_zcoordinate.append(rnew_ycoordinate[a]*math.sin(math.radians(beta)) + rnew_zcoordinate[a]*math.cos(math.radians(beta)))

	for a in range(0,natom):
		if p2[a] == 1 :        #a gives index of the unique nitrogen which is connected to one carbon atom. We are rotating the molecule to bring this atom in the yz plane
			print a+1,atom_type[a], p2[a]
			x3=r2new_xcoordinate[a]
			y3=r2new_ycoordinate[a]
			z3=r2new_zcoordinate[a]
			ceta = fsolve(rotation_y,(0))
			atom3=a

	z_atom3 = r2new_zcoordinate[atom3]*math.cos(math.radians(ceta)) - r2new_xcoordinate[atom3]*math.sin(math.radians(ceta))   #z coordinate after rotation along y axis

	if z_atom3 < 0:   #if rotation has taken place such that the third nitrogen sits in negative yz plane
		ceta=ceta+180   #rotate in the opposite direction such that the third nitrogen is towards the positive z direction in yz plane
	for a in range(0,natom):
		r3new_xcoordinate.append(r2new_xcoordinate[a]*math.cos(math.radians(ceta)) + r2new_zcoordinate[a]*math.sin(math.radians(ceta)))    #new and final coordinates of the molecule after reorientation
		r3new_ycoordinate.append(r2new_ycoordinate[a]*1 )
		r3new_zcoordinate.append(r2new_zcoordinate[a]*math.cos(math.radians(ceta)) - r2new_xcoordinate[a]*math.sin(math.radians(ceta)))

	for a in range(0,natom):
		dist_origin.append(dist[atom1][a])    #matrix which gives distance of all atoms from the origin
	
	return (natom, atom_type, r3new_xcoordinate, r3new_ycoordinate, r3new_zcoordinate, dist_origin) 

natom_ref, atom_type_ref, xcoordinate_ref, ycoordinate_ref, zcoordinate_ref, dist_ref = orientation(file_ref)    #reference file
natom_2, atom_type_2, xcoordinate_2, ycoordinate_2, zcoordinate_2, dist_2 = orientation(file_2)                  #file to be renumbered. Should have larger or same number of atoms than reference file

#Now the renumbering begins :)

X = [[0 for i in range(natom_2)] for j in range(natom_2)]    #2D matrices of dimension natom_2 X natom_2   given natom_2 > natom_ref
Y = [[0 for i in range(natom_2)] for j in range(natom_2)]
Z = [[0 for i in range(natom_2)] for j in range(natom_2)]

new_order = []
prev_index=[]
xnew = []
ynew = []
znew = []

for a in range(natom_ref,natom_2):
	atom_type_ref.append('X')             #putting 'X' in atom type of reference file for the extra atoms in file2
	xcoordinate_ref.append(0.0)           #putting 0.0 in x coordinate of reference file for the extra atoms in file2
	ycoordinate_ref.append(0.0)           #putting 0.0 in y coordinate of reference file for the extra atoms in file2
	zcoordinate_ref.append(0.0)           #putting 0.0 in z coordinate of reference file for the extra atoms in file2
	dist_ref.append('nan')                #putting no value in distance from origin of reference file for the extra atoms in file2


for a in range(0,natom_2):
	for b in range(0,natom_2):       #how distant each point in reference file is from each point in file2
		X[a][b]	= xcoordinate_ref[a] - xcoordinate_2[b]
		Y[a][b]	= ycoordinate_ref[a] - ycoordinate_2[b]
		Z[a][b]	= zcoordinate_ref[a] - zcoordinate_2[b]

for a in range(0,natom_2):
	for b in range(0,natom_2):
		if atom_type_ref[a] != 'X' :
			if (abs(X[a][b]) < 0.19) and (abs(Y[a][b]) < 0.19) and (abs(Z[a][b]) < 0.19):    #if any atom in reference file has the same position in file2, give the same number to the atom
				new_order.append(atom_type_2[b])
				xnew.append(xcoordinate_2[b])
				ynew.append(ycoordinate_2[b])
				znew.append(zcoordinate_2[b])	
				prev_index.append(b)


for a in range(0,natom_2):
	for b in range(0,natom_2):
		if a not in prev_index:          #in case file2 has more atoms than reference file, number the extra atoms in a random way after the matching atoms are renumbered 
			new_order.append(atom_type_2[a])
			xnew.append(xcoordinate_2[a])
			ynew.append(ycoordinate_2[a])
			znew.append(zcoordinate_2[a])	
			prev_index.append(a)
#		if not atom_type_ref[a] != 'X' :
#			new_order.append(atom_type_2[b])
#			xnew.append(xcoordinate_2[b])
#			ynew.append(ycoordinate_2[b])
#			znew.append(zcoordinate_2[b])	
		
for a in range(0,natom_2):
	print new_order[a],xnew[a],ynew[a],znew[a] #atom_type_ref[a]
	#print prev_index[a]

