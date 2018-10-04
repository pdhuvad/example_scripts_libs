#!/usr/bin/python

#####################################################################
#####################################################################
# Script to calculate total polarization using born effective charge#
#####################################################################
#####################################################################


#####################################################################
###########                INSTRUCTIONS                 #############
# Files required
# 1) OUTCAR file with Born-effective charges
# 2) CONTCAR   : low-symmetry CONTCAR file
# 3)** CONTCAR_high  : high-symmetry CONTCAR file  

#  ** if high symmetry positions geneterated by this script are 
#     incomensurate then  make appropriate changes ( uncomment #### lines)



###################         OPTIONS                 #################
## if DIAG = yes ; use only diagonal elements of Born-effective
## charges to calcculate polarization   

#####################################################################


import numpy as np

#DIAG = 'yes'
DIAG = 'no'

np.set_printoptions(precision=4,formatter={'float': '{: 0.9f}'.format})

def high_sym_pos(pos):
	if (pos >= -.125 and pos <= .125):
		pos_hisym = 0.0
	elif (pos >= .125 and pos <= .375):
		pos_hisym = 0.25
	elif (pos >= .375 and pos <= .625):
		pos_hisym = 0.5
	elif (pos >= .625 and pos <= .875):
		pos_hisym = 0.75
	return pos_hisym


def convert_direct(tmp):
	for i in range(len(tmp)):
	    if tmp[i][0] > .9375:
	        tmp[i][0]= tmp[i][0] -1
	    for j in [1,2]:
	        if tmp[i][j] > .875:
	            tmp[i][j] = tmp[i][j] -1
	            
	return tmp


file_born = open('OUTCAR',"r")
file_cls = open('CONTCAR','r')

position_ls=[]
for line in file_cls:
	ln = line.split()
	position_ls.append(ln)
prim_mat=position_ls[2:5]
prim_mat=[[float(k) for k in l] for l in prim_mat]
N = int(sum([float(k) for k in position_ls[6]]))
tmp_position = position_ls
dum = N + 8
position_ls = position_ls[:dum]
position_ls = position_ls[8:]
position_ls = [[float(i) for i in k] for k in position_ls]
position_ls = np.array(convert_direct(position_ls))
print "####################### low symmetry positions #########################"
print position_ls




a=[]
for i in range(8):
    b= .125*i
    a.append(b)
print "a",a
intervals=[]
for i in a:
    intervals.append([i, i-.0625 , i+.0625, ])
print "intervals", intervals


position_hs=np.copy(position_ls)
print "position_hs",position_hs
position_hs=convert_direct(position_hs)
for k in range(len(position_hs)):
        print k
        i =position_hs[k][0]
        for j in range(len(intervals)):
                if (i >= intervals[j][1]) and ( i <= intervals[j][2]) :
                        dum = intervals[j][0]
                        position_hs[k][0] = np.copy(dum)
                else:
                        pass
        position_hs[k][1] = high_sym_pos(position_hs[k][1])
        position_hs[k][2] = high_sym_pos(position_hs[k][2])
print "position_hs",position_hs

print "########### script generated high symmetry positions ###################"
#print position_hs

for i in range(len(position_hs)):
	print position_ls[i][0] , "--->" , position_hs[i][0] , "\t" , position_ls[i][1] , "--->" , position_hs[i][1] , "\t" , position_ls[i][2] , "--->" , position_hs[i][2]
####file_chs = open('CONTCAR_high','r')
####position_hs=[]
####for line in file_chs:
####	ln = line.split()
####	position_hs.append(ln)
####dum = N + 8
####position_hs = position_hs[:dum]
####position_hs = position_hs[8:]
####position_hs = [[float(i) for i in k] for k in position_hs]
###
####position_hs = np.array(convert_direct(position_hs))
atom_disp = position_hs - position_ls
for i in range(N):
	atom_disp[i] = np.matmul(prim_mat,atom_disp[i])

born_matrix=[]
for line in file_born:
	if '  volume of cell : ' in line:
		ln = line.split()
		volume = ln[4]
		volume = float(volume)
		
	if 'BORN EFFECTIVE CHARGES ' in line:
		file_born.next()
		for i in range(N):
			born_array=[]
			#print "i is",i
			for j in range(4):
				ln = file_born.next()
				if 'ion' in ln:
					pass
				else:
					row=[]
					row=ln.split()
					row = row[1:]
					row = [float(k) for k in row ]
					born_array.append(row)
			born_matrix.append(born_array)

bm_diag = np.copy(born_matrix)
bm_diag
for i in range(len(bm_diag)):
    for j in range(3):
        for k in range(3):
            if j ==k:
                pass
            else:
                bm_diag[i][j][k] = 0

if DIAG == 'yes':
	born_matrix = np.copy(bm_diag)
else:
	pass

atomic_dipol = []
for i in range(N):
    row = []
    row = np.matmul(born_matrix[i],atom_disp[i])
    atomic_dipol.append(row)
    
atomic_dipol = np.array(atomic_dipol)
#print "volume", volume
total_dipol = [0,0,0]
for i in atomic_dipol:
    total_dipol = np.add(total_dipol, i)
    
total_polarization = 1602 * total_dipol/volume
print "total polarization in uC/cm^2" , total_polarization
                
