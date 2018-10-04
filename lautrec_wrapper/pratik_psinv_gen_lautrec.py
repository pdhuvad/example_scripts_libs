#!/usr/bin/python
import numpy as np
import con_rel  # module to generate coo_start and check force convergance
import os
from  subprocess import Popen
import fileinput
import shutil




#initialization tags
lautrec = 'echo "-n 256 /global/homes/p/pdhuvad/LAUTREC/compile/lautrecPAW-1.0.0_x86_64.x xx"'
delta = 0.01   #dispacement given to each dof 
wannierization = 1    # 1=on , 0= off
D_list = str(.015 *8.6467749549)
constd = '.t.'
init_wvfn = 'true'
init_loop = 'true'
zz=[0]*31
for i in range(30):
	print zz[i]
zz[0]=0
zz[1] =  0.007474
zz[2] =  0.692917
zz[3] =  0.348009
zz[4] =  1.034601
zz[5] = -0.010648
zz[6] =  0.388442
zz[7] =  0.273152
zz[8] =  0.674289
zz[9] =  0.959013
zz[10] =  1.075254
zz[11] =  0.030481
zz[12] = -0.043449
zz[13] = -0.001423
zz[14] =  0.000427
zz[15] = -0.039624
zz[16] =  0.191370
zz[17] =  0.263431
zz[18] = -0.001534
zz[19] =  0.202125
zz[20] =  0.267806
zz[21] =  0.509471
zz[22] =  0.494166
zz[23] =  0.000308
zz[24] =  0.004852
zz[25] =  0.082066
zz[26] =  0.285073
zz[27] =  0.786312
zz[28] = -0.080558
zz[29] =  0.278352
zz[30] =  0.787044
print zz
zz=np.array(zz)
ndof=int(len(zz))
print "number of degrees of freedom", ndof
def get_rotations(zz):


	xx=[0]*11; yy=[0]*11
	
	xx[1] = zz[11]+0.5
	xx[2] = zz[12]+0.5
	xx[3] = zz[13]+0.5
	xx[4] = zz[14]+0.5
	xx[5] = zz[15]+0.5
	xx[6] = zz[16]+0.5
	xx[7] = zz[17]+0.5
	xx[8] = zz[18]+0.5
	xx[9] = zz[19]+0.5
	xx[10] = zz[20]+0.5
	yy[1] = 0.5-zz[21]
	yy[2] = 0.5-zz[22]
	yy[3] = 0.5-zz[23]
	yy[4] = 0.5-zz[24]
	yy[5] = 0.5-zz[25]
	yy[6] = 0.5-zz[26]
	yy[7] = 1.5-zz[27]
	yy[8] = 0.5-zz[28]
	yy[9] = 0.5-zz[29]
	yy[10] = 1.5-zz[30]
	xcor = xx
	ycor = yy
	print xx, yy
	return (xx, yy)

xcord,ycord = get_rotations(zz)
print "xcord, ycord",xcord,ycord
def get_xred(xx,yy,zz):
	
	xred=[]
###gaur###
	xred.append( [zz[1] ,  zz[11] ,  zz[21]  ] )  
	xred.append( [zz[1] ,  xx[1]  ,  yy[1]   ] )
	xred.append( [zz[2] ,  zz[12] ,  zz[22]  ] )
	xred.append( [zz[2] ,  xx[2]  ,  yy[2]   ] )
	
	xred.append( [zz[3] ,  zz[13] ,  zz[23]  ] )
	xred.append( [zz[3] ,  xx[3]  ,  yy[3]   ])
	xred.append( [zz[4] ,  zz[14] ,  zz[24]  ] )
	xred.append( [zz[4] ,  xx[4]  ,  yy[4]   ])
	
	xred.append( [zz[5] ,  zz[15] ,  zz[25]  ] )
	xred.append( [zz[5] ,  xx[5]  ,  yy[5]   ])
	xred.append( [zz[6] ,  zz[16] ,  zz[26]  ] )
	xred.append( [zz[6] ,  xx[6]  ,  yy[6]  ])
	
	xred.append( [zz[7] ,  zz[17] ,  zz[27]  ] )
	xred.append( [zz[7] ,  xx[7]  ,  yy[7]   ])
	xred.append( [zz[8] ,  zz[18] ,  zz[28]  ] )
	xred.append( [zz[8] ,  xx[8]  ,  yy[8]   ])
	
	xred.append( [zz[9] ,  zz[19] ,  zz[29]  ] ) 
	xred.append( [zz[9] ,  xx[9]  ,  yy[9]   ] )
	xred.append( [zz[10],  zz[20] ,  zz[30]  ] ) 
	xred.append( [zz[10],  xx[10] ,  yy[10]  ] )
###nitai###	
	xred=np.array(xred)
	return xred


def get_force(nat,ndof):   #This function uses "rotation" and "out_forces" files to get mapped forces for degree of freedom.
			   # mapping is done beteween rotations and out_forces force.
			   # Here the forces are averaged for given degree of freedom for corrospoinding rotation index. 
	
	file_force=open("out_forces","r")
	for line in file_force:
	
		if line == '  *** Coordinates and forces ***\n':
			force_lautrec=[]
			for _ in range(nat):
				ln= file_force.next()
				row=[]
				row=ln.split()
				lst=[]
				for j in [5,6,7]:   # last three columes are forces
					const=float(row[j])  #here put check for float, if not exit and show error
					lst.append(const)
				force_lautrec.append(lst)
	
	file_force.close()
	force_lautrec=np.array(force_lautrec)
	print "force_lautrec", force_lautrec
	file_data=open("rotation","r")
	i=0
	mapped_force={}
	for line in file_data:
		ln=str.split(line)
		print "ln", ln
		for j in range(3):
			mapped_force.setdefault(str(ln[j]),[]).append(force_lautrec[i][j])  #mapping
		i += 1
	file_data.close()
	zz_force=[]
	for i in range(1,ndof):
		num="".join(["zz[",str(i),"]"])
		num1 =0
		for j in range(len(mapped_force[num])):
			num1=num1 + float(mapped_force[num][j])
		zz_force.append(num1)
		print num, num1, mapped_force[num]
	zz_force=np.array(zz_force)
	print zz_force
	return zz_force
def get_alat():
	inp10_file=open("xx_inp10","r")
	i = 0
	for line in inp10_file:
		ln=line.split()
		print ln
		alat=0.
		if i==0:
			alat=float(ln[0])
			print "alat in loop", alat
			break
	return alat

os.system("tail -n 6 pratik_psinv_gen_lautrec >tmpfile1")
os.system("sed 's/###123//g' tmpfile1>tmpfile2")
os.system("bash tmpfile2")

def copyLargeFile(src, dest, buffer_size=16000):
	with open(src, 'rb') as fsrc:
		with open(dest, 'wb') as fdest:
			shutil.copyfileobj(fsrc, fdest, buffer_size)

############################ Initialization of wavefunction#####################
file_log= open("log_run", "w+")

## create inp10_98 and tail_10 files from xx_inp10   , this is unnecessacery , I dont know why I did it#####
file_tail_10 = open("tail_10",'w+',0)
file_inp10 = open("xx_inp10",'r')
file_inp10_98= open("inp10_98","w+",0)
k = 0
for line in file_inp10:
	k = k + 1
	ln=line
	if k < 5 :
		file_inp10_98.write(ln)
	if k > 4:
		file_tail_10.write(ln)
file_inp10.close()
file_inp10_98.close()
file_tail_10.close()
###########################################################################

#os.system("cat inp10_98 tail_10 > xx_inp10")
#os.system("./tiling_4kpar.x xx 16")


############### create a series of inp20_98 files #################
fname=str('inp20_98_DOF' + str(0))
file_inp20_98=[0]*ndof
savetxt_file_inp20_98=[0]*ndof
int_para = get_xred(xcord,ycord,zz)
print "int_para", int_para
xx0,yy0=get_rotations(zz)
original_int_para=get_xred(xx0,yy0,zz)
natoms=len(get_xred(xx0,yy0,zz))
#### column of 1s which is imove
new_col=np.ones((natoms,1),dtype='int')
print new_col

####################

forcebefore=[]; forceafter=[]
file_inp20_98[0]=open(fname, "w")
np.savetxt(file_inp20_98[0],np.append(original_int_para,new_col,1),fmt='%1.8f \t  %1.8f \t %1.8f \t %d',delimiter='\t')
file_inp20_98[0].close()
int_para=[0]*ndof
for i in range(1,ndof):
	new_zz = np.copy(zz)
	new_zz[i]=zz[i]+delta
	xcord,ycord = get_rotations(new_zz)
	int_para[i]=get_xred(xcord,ycord,new_zz)
	print "lengths" ,len(int_para[i]), len(new_col)
	print "int_para", i 
	print int_para[i] 
	#print original_int_para
	print "savetxt",np.append(int_para[i],new_col,1)
	fname=str('inp20_98_DOF' + str(i))
	file_inp20_98[i]=open(fname, "w")
	np.savetxt(file_inp20_98[i],np.append(int_para[i],new_col,1),fmt='%1.8f \t  %1.8f \t %1.8f \t %d',delimiter='\t')
	file_inp20_98[i].close()
	if init_wvfn == 'true' and i == 1:
		print "Here"

		copyLargeFile("inp20_98_DOF0","inp20_98")
		file_log.write('Initializing wavefunction... \n')
		file_inp20_00 = open("xx_inp20_00",'w+',0)
		file_inp20_01 = open("xx_inp20_01",'w+',0)
		for line in fileinput.input(['head_00','inp20_98']):
		        file_inp20_00.write(line)
		file_inp20_00.close()
		for line in fileinput.input(["head_91","inp20_98"]):
		        file_inp20_01.write(line)
		file_inp20_01.close()
	
		os.system("./inpgen 11 12 .f. 0.0 ")
	
		print "00"
		file_log.write("00 started")
		lautrec_00 = ' '.join([lautrec,"00", ">out_00"])
		file_log.write("LAUTREC 00 > out_00  if any stderr ...\n")
		os.system(lautrec_00)
	
		print "01"
		file_log.write("01 started")
		lautrec_01 = ' '.join([lautrec,"01", ">out_01"])
		file_log.write("LAUTREC 01 > out_01  if any stderr ...\n")
		os.system(lautrec_01)
	
		print "11"
		lautrec_11 = ' '.join([lautrec,"11", ">out_11"])
		file_log.write("LAUTREC 11 > out_11  if any stderr ...\n")
		os.system(lautrec_11)
	
		print "12"
		lautrec_12 = ' '.join([lautrec,"12", ">out_12"])
		file_log.write("LAUTREC 12 > out_12  if any stderr ...\n")
		os.system(lautrec_12)
	
	#Initialize the loop
	
	if init_loop == 'true' and i == 1:
		file_log.write('Initializing the loop..\n')
		inp_dfield=' '.join(["./inpgen 13 14", constd, D_list])
		os.system(inp_dfield)
	
		print "13"
		file_log.write('init_loop 13 started \n')
		lautrec_13 = ' '.join([lautrec,"13", ">out_13"])
		file_log.write("LAUTREC 13 > out_13  if any stderr ...\n")
		os.system(lautrec_13)
	
		print "14"
		file_log.write('init_loop 14 started \n')
		lautrec_14 = ' '.join([lautrec,"14", ">out_14"])
		file_log.write("LAUTREC 14 > out_14  if any stderr ...\n")
		os.system(lautrec_14)
		os.system("cp out_14 out_forces")
		os.system("cp out_forces_00 out_forces")
		nat = 20
		forcebefore_array=[]
		forcebefore_array=get_force(nat,ndof)
		forcebefore_array=np.array(forcebefore_array)
		for i in range(1,ndof):
			forcebefore.append(forcebefore_array)
		forcebefore=np.array(forcebefore)
		file_log.write(" \n\n forcebefore matrix is generated \n\n")
	file_log.write(" Begin total energy calculations for perturbed DOF\n")
	copy_inp20_98_dof = ' '.join(["cp", fname,"inp20_98"])
	os.system(copy_inp20_98_dof)

	inp_dfield_dof=' '.join(["./inpgen 01 02", constd, D_list])
	os.system(inp_dfield_dof)
	
	file_log.write('TOTAL ENERGY Calculation for DOF No: {0:2d}  \n'.format(i))
	lautrec_dof = ' '.join([lautrec,"01 >",''.join(["out-dof","-01-",str(i)]) ])
	file_log.write('01  if any stderr ...\n')
	os.system(lautrec_dof)
	
	lautrec_dof = ' '.join([lautrec,"02 >",''.join(["out-dof","-02-",str(i)]) ])
	file_log.write('02  if any stderr ...\n')
	os.system(lautrec_dof)
	
	cp_out_string=''.join([''.join(["cp out-dof","-02-",str(i)]), " out_forces"])
	os.system(cp_out_string)
	os.system("cp out_forces_00 out_forces")
	forceafter.append(get_force(nat,ndof))
forceafter = np.array(forceafter)
file_log.write(" \n\n forceafter matrix is generated \n\n")
print "forceafter" , forceafter 
print "forcebefore", forcebefore
print "forcebefore_array", forcebefore_array
fb_file=open("forcebefore","w")
fa_file=open("forceafter","w")
np.savetxt(fb_file, forcebefore,fmt="%1.8f")
np.savetxt(fa_file, forceafter,fmt="%1.8f")
fb_file.close()
fa_file.close()
file_log.close()
forcebefore=np.random.randn(30,30)
forceafter=np.random.randn(30,30) * 2

diff_force=forceafter - forcebefore
print get_alat(), delta
perturbation = -1 * delta * get_alat()
file_log.write(" \n perturbation in cartezian coordinates, delta , alattice %1.6f %1.6f %1.5f \n\n"%(perturbation, delta, get_alat()))
force_constant = diff_force/perturbation
eigenvalues_force=np.linalg.eig(force_constant)
print "eigenvalues", eigenvalues_force
psinv_mat=np.linalg.pinv(force_constant)
print "psinv_mat", psinv_mat
psinv_file=open("psinv.dat","w")
np.savetxt(psinv_file,psinv_mat,fmt="%1.10f")
psinv_file.close()

print "number of degrees of freedom", ndof





###123#!/bin/sh
###123awk '/###gaur/{flag=1;next}/###nitai/{flag=0} flag{print}' pratik_psinv_gen_lautrec >tmp1
###123cat tmp1 | awk '{gsub(/[[:blank:]]/, "");gsub("xred.append\\(\\[","");gsub("\\]\\)","");gsub(",","\t");gsub("\\n","");print}' >tmp2
###123awk 'NF >0' tmp2 > tmp3
###123sed '$d' tmp3 |sed '$d' | sed '$d' >rotation ; rm tmp*
