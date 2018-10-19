#!/usr/bin/python
##########################################################################
#################### Python wraper for LAUTREC ###########################
# __auther__ : Pratik Dhuvad
# __version__ : 0.2
# __build__ : July 2014
# __contact__ : tud12582@temple.edu
# 
#
#   Changes:  
#       a)  automated initialization of starting coordinates
#       b) force convergence check for each BFGS loop
#       
#########################################################################
import con_rel  # module to generate coo_start and check force convergance
import numpy as np
import os
from  subprocess import Popen
import fileinput
#control parameters
elastic_constant = .52/60
nat = 40  # no of atoms 

# TODO: can be parsed from the coo_start

start_from_scratch = 'yes'

#BFGS loop control parameters
init_coords = 'true'
init_wvfn = 'true'
init_loop =  'true'
init_hessian = 'old'

## TODO:  parse these parameters from argument

max_inner_iterations =35
max_outer_iterations =7
pulay = '0.000'
D_list = np.array([        0.055  ,   0.060  ,   0.065  ,   0.070])
D_lautrec_coords = 10.379225950000*10.379225950000/(4*3.14159265358979323846)
#d_list = [str(x) for x in d_list]
#D_list = str( 0.000 * 10.37922695898 * 10.37922695898/(4*3.14159265358979323846))
constd = '.t.'

#LAUTREC executable 
#lautrec = "srun -n 256 /global/homes/p/pdhuvad/LAUTREC/compile/lautrecPAW-1.0.0_x86_64.x xx"
#lautrec = "srun -n 256 /global/homes/p/pdhuvad/LAUTREC/old_lautrec/lautrecPAW-1.0.0_x86_64.x xx"
#lautrec = "srun -n 256 /global/u2/w/whw1985/Code/LAUTREC_mohan/lautrecPAW-1.0.0_i386.x xx"
lautrec = "mpirun -np 64  /home/tud12582/LAUTREC/Lautrec/pbesol_pratik/lautrecPAW-1.0.0_x86_64.x  xx"
#lautrec = "mpirun -np 64 --map-by socket lautrecPAW-1.0.0_x86_64.x xx " 
#lautrec = "srun -n 256 /global/homes/p/pdhuvad/LAUTREC/hongwei/LAUTREC_mohan/lautrecPAW-1.0.0_i386.x xx"
lautrec_test = 'echo "-n 256 /global/homes/p/pdhuvad/LAUTREC/compile/lautrecPAW-1.0.0_x86_64.x xx"'
with open ("log_run",'w+',0) as file_log:
	with open ("int_para",'w+',0) as file_int_para:
			# D_list loop_starts here
		for k in D_list:
			d_value = k*D_lautrec_coords
                        file_log.write("#################################################################################\n")
                        file_log.write("####################### D_field {0:<8.4f}  ###################################\n".format(k))
                        file_log.write("#################################################################################\n")
		#Outer loop
		
			for j in range(1,max_outer_iterations):
				if start_from_scratch == 'no':
					file_log.write("*********************************************************************************\n")
					file_log.write("start_from_scratch = no , using old hessian and wavefunction for first outer step\n")
					file_log.write("*********************************************************************************\n")
				file_log.write("*********************************************************************************\n")
				file_log.write("*********************** Outer ITERATION {0:2d} ***********************************\n".format(j))
				file_log.write("*********************************************************************************\n")
				
				file_int_para.write("*********************************************************************************\n")
				file_int_para.write("*********************** Outer ITERATION {0:2d} ***********************************\n".format(j))
				file_int_para.write("*********************************************************************************\n")
				
				#initializing the coordinates
				if init_coords == 'true' and start_from_scratch =='yes':
					file_log.write ('Initializing coordinates... \n')
					file_coo_start = open("coo_start",'r')
					file_inp10_98 = open("inp10_98",'w+',0)
					file_inp20_98 = open("inp20_98",'w+',0)
					file_tail_10 = open("tail_10",'w+',0)
					file_inp10 = open("xx_inp10",'r')
					#for i in range(4):
					#	file_inp10_98.write(file_coo_start.next())
					kl = 0
					for line in file_coo_start:
						kl = kl + 1
						ln=line
						if kl < 5 :
							file_inp10_98.write(ln)
						if kl == 5:
							pass
						if kl > 5:
							file_inp20_98.write(ln)
					kl = 0
					for line in file_inp10:
						kl = kl + 1
						ln=line
						if kl < 5 :
							pass
						if kl > 4:
							file_tail_10.write(ln)
					file_tail_10.close()
					file_inp10_98.close()	
					file_inp10_98.close()	
					file_inp20_98.close()	
					file_coo_start.close()
					os.system("cat inp10_98 tail_10 > xx_inp10")
					os.system("./tiling_p4mm.x xx 16")       
				# Initialize the wavefunction
				if init_wvfn == 'true' and start_from_scratch =='yes':
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
				
				if init_loop == 'true' and start_from_scratch =='yes' :
					file_log.write('Initializing the loop..\n')
					inp_dfield=' '.join(["./inpgen 13 14", constd, str(d_value)])
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
				
				
				#Initialize the Hessian
				if init_hessian == 'new' and start_from_scratch =='yes' and k == D_list[0]:
					inp_1_usefmat=' '.join(["echo 1",str(d_value), pulay,"| ./usefmat-non-sym.x >", ''.join(["log-",str(k),"-0"])]) #,''.join(">log-",d_list,"-1"))
					os.system(inp_1_usefmat)
					os.system("cp hessian.dat hessian.dat-1")
				
				if init_hessian == 'old':
					os.system("cp hessian.dat_old hessian.dat")
					inp_0_usefmat=' '.join(["echo 0",str(d_value), pulay,"| ./usefmat-non-sym.x >",''.join(["log-",str(k),"-1"])])
					os.system(inp_0_usefmat)
					os.system("cp hessian.dat hessian.dat-0")
				
				
				#BFGS loop starts here
				start_from_scratch ='yes'
				for i in range(11,max_inner_iterations):
					
					file_log.write('BFGS/Inner loop..\n')
					inp_dfield_bfgs=' '.join(["./inpgen 01 02", constd, str(d_value)])
					os.system(inp_dfield_bfgs)
				
					file_log.write('BFGS iteration No: {0:2d}  \n'.format(i))
					lautrec_bfgs = ' '.join([lautrec,"01 >",''.join(["out-",str(k),"-01-",str(i)]) ])
					file_log.write('01  if any stderr ...\n')
					os.system(lautrec_bfgs)
				
					lautrec_bfgs = ' '.join([lautrec,"02 >",''.join(["out-",str(k),"-02-",str(i)]) ])
					file_log.write('02  if any stderr ...\n')
					os.system(lautrec_bfgs)
					
					cp_out_string=''.join([''.join(["cp out-",str(k),"-02-",str(i)]), " out_forces"])
					os.system(cp_out_string)
        	                        con_rel.make_coo_start()
        	                        file_log.write('coo_start has been updated\n')
                                	os.system("python symmetrize_coo_start.py")
                                	file_log.write("coo_start has been symmetrized")
					inp_bfgs_usefmat=' '.join(["echo",str(i),str(d_value), pulay,"| ./usefmat-non-sym.x ",''.join([">log-",str(k),"-",str(i)])])
					os.system(inp_bfgs_usefmat)
					cp_hess_string=''.join(["cp hessian.dat hessian.dat-",str(k),"-",str(i)])
					os.system(cp_hess_string)
					os.system("cp hessian.dat hessian.dat_old")
					cp_vprev_string=''.join(["cp vprev.dat vprev.dat",str(k),"-",str(i)])
					os.system(cp_vprev_string)
					os.system("cp out_forces OUT")
				#copy primitive vectors, force and stress in a file
					file_int_para.write("*********************************************************************************\n")
					file_int_para.write("*********************** BFGS ITERATION {0:2d} ***********************************\n".format(i))
					file_int_para.write("*********************************************************************************\n")
				        f =  open('out_forces','r')
				        prim_mat=[]
				        print "hari bol"
				        # folowing for loop reads lattice parameters, internal parameters and stress
				        for line in f:
				                if line ==' Real primitive vectors (atomic units): \n':
							file_int_para.write(' Real primitive vectors (atomic units): \n')
				                        for _ in range(3):
				                                ln= f.next()
								file_int_para.write(ln)
				                if line == '  *** Coordinates and forces ***\n':
							file_int_para.write('  *** Coordinates and forces ***\n')
				                        for _ in range(nat):
				                                ln= f.next()
								file_int_para.write(ln)
				                if line == ' *** TOTAL stress (before) ***\n':
							file_int_para.write(' *** TOTAL stress (before) ***\n')
				                        for _ in range(3):
				                                ln= f.next()
								file_int_para.write(ln)
				#                if line == '  *** Coordinates and forces ***\n':
				#			list_force_ln=[]
				#                        for _ in range(nat):
				#                                ln= f.next()
				#				list_force_ln.append(ln)
				#                if line == ' *** TOTAL stress (before) ***\n':
				#			list_stress_ln=[]
				#                        for _ in range(3):
				#                                ln= f.next()
				#				list_stress_ln.append(ln)
				#	file_int_para.write(' *** TOTAL stress (before) ***\n')
				#	file_int_para.write(list_stress_ln)
				#	file_int_para.write('  *** Coordinates and forces ***\n')
				#	file_int_para.write(list_force_ln)
					file_int_para.flush()
				
				#check for force convergance... 
					const=con_rel.check_force()
					if const == True:
				        	file_log.write( "the forces are not converged, the BFGS will continue\n")
					else:
				        	file_log.write("force is converged, end BFGS\n")
				        	break
					
					
				con_rel.make_coo_start()	
				os.system("python symmetrize_coo_start.py")
				file_log.write("coo_start has been symmetrized")
				os.system("cp hessian.dat hessian.dat_old")
				#check for stress convergance
				con_str=con_rel.check_stress()
				init_hessian = 'old'
				if con_str == True:
			        	file_log.write( "stress is  not converged, the outer loop  will continue\n")
				        #os.system("rm grid* coord.* Timing* wvfn* fort* xx_inp20* stdout* vhart* ")
				else:
			        	file_log.write("stress is converged, calculations ENDS\n")
			        	break
		
