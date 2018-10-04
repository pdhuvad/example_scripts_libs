#!/usr/bin/python
import numpy as np
import sys
elastic_constant = .65/52 #elastic constant of CaTiO3 P4bm
strain_step = 0.008
min_force = 0.000001
min_strain = 0.000001
nat = 20 #number of atoms
global strs_x
def make_coo_start():
	#parameters
	f =  open('out_forces','r')
	prim_mat=[]
	print "hari bol"
	# folowing for loop reads lattice parameters, internal parameters and stress
	for line in f:
		if line ==' Real primitive vectors (atomic units): \n':
			for _ in range(3):
				ln= f.next()
				row=[]
				row=ln.split()
				for j in [0,1,2]:
					row[j]=float(row[j])  #here put check for float, if not exit and show error
				prim_mat.append(row)
		#	prim_mat=np.array(prim_mat)
			print prim_mat
			#prim_mat=np.matrix(np.array(prim_mat))
			prim_mat=np.array(prim_mat)
		if line == '  *** Coordinates and forces ***\n':
			cart_cord=[]
			for _ in range(nat):
				ln= f.next()
				row=[]
				row=ln.split()
				lst=[]
				for j in [2,3,4]:
					const=float(row[j])  #here put check for float, if not exit and show error
					lst.append(const)
				cart_cord.append(lst)
		if line == ' *** TOTAL stress (before) ***\n':
			stress=[]
			for _ in range(3):
				ln= f.next()
				row=[]
				row=ln.split()
				lst=[]
				for j in [0,1,2]:
					const=float(row[j])  #here put check for float, if not exit and show error
					lst.append(const)
				stress.append(lst)
	stress=np.array(stress)
	cart_cord=np.array(cart_cord)
	print prim_mat
	print cart_cord
	print 'stress: \n  ', stress
	print np.linalg.inv(prim_mat)
	
	frac_cord=[]
	print '\n\n\n'
	for i in range(0,len(cart_cord)):
		row=[]
		row=np.dot(np.linalg.inv(np.transpose(prim_mat)),cart_cord[i])
		frac_cord.append(row)
		#print row
	frac_cord=np.array(frac_cord)
	print frac_cord
	
	#change of cell parameters
	lattice_x = prim_mat[0][0]
	stress_x= stress[0][0]
	strain = stress_x/elastic_constant
	print "strain=" , strain, "stress_x", stress_x, "lattice_x", lattice_x
	if abs(strain) > strain_step:
		strain = strain / abs(strain) * strain_step
	lattice_x = lattice_x *( 1 + strain)
	prim_mat[0][0] = lattice_x
	new_cart_cord=[]
	for i in range(0,len(cart_cord)):
		row=[]
		row=np.dot(np.transpose(prim_mat),frac_cord[i])
		new_cart_cord.append(row)	
	
	### print coo_start file
	new_col=np.ones((nat,1),dtype='int')
	coo_start_new_cart_cord=np.append(new_cart_cord,new_col,1)
	file_coo_start=open("coo_start", 'w+')
	file_coo_start.write('    1.0 \n')
	np.savetxt(file_coo_start,prim_mat,fmt='%.8f',delimiter='\t')
	file_coo_start.write('----------------------------------------------------------------\n')
	np.savetxt(file_coo_start,coo_start_new_cart_cord,fmt='%2.10f          %2.10f        %2.10f\t %1d',delimiter='\t')


#'Check force' module for inner BFGS loop
# access this function from bash as: python -c 'import con_rel; print con_rel.check_force()'
# if _ = True then continue
# if _ = None converged, exit BFGS 
def check_force():
	f =  open('out_forces','r')
	for line in f:
		if line == '  *** Coordinates and forces ***\n':
			forces=[]
			for _ in range(nat):
				ln= f.next()
				row=[]
				row=ln.split()
				lst=[]
				for j in [5,6,7]:
					const=float(row[j])  #here put check for float, if not exit and show error
					lst.append(const)
				forces.append(lst)
	forces_array=np.array(forces)
	#print forces
	for list in forces_array:
		for element in list:
			if abs(element) > min_force:
				return True
# this function checks for the stress convergance

def check_stress():
	f = open('out_forces','r')
	for line in f:
		if line == ' *** TOTAL stress (before) ***\n':
			ln=[]
			for _ in range(1):
	                	ln=f.next()
				strs_row=ln.split()
				strs_x = float(strs_row[0])
				break
	
	strs_x= strs_x/elastic_constant
	if abs(strs_x) > min_strain:
		return True
