#!/usr/bin/python
import numpy as np

##### TO Read OUTCAR file in the same direcotyr##########
f= open("OUTCAR", "r")
for line in f:
        if "p[ion]=(" in line:
                ln= line.split()
                #print ln
                ionic_contribution=[]
                for i in [4,5,6]:
                        ionic_contribution.append(float(ln[i]))
######################### NOTE ###########################
## for spin polarized calculations , only total electrionc
## contribution is considered
        if "p[elc]=(" in line:   
                ln1= line.split()
                #print ln1
                ele_contribution=[]
                for i in [5,6,7]:
                        ele_contribution.append(float(ln1[i]))
        if "volume of cell :" in line:
                ln1= line.split()
        #       #print ln1
                omega=float(ln1[4])

        if "direct lattice vectors" in line:
                prim_mat=[]
                for i in range(3):
                        ln=f.next()
                        ln = ln.split()
                        #print ln[0]
                        row=[]
                        for i in range(3):
                                row.append(float(ln[i]))
                                #print ln[i]
                        prim_mat.append(row)
f.close()
ionic_contribution = np.array(ionic_contribution)
ele_contribution = np.array(ele_contribution)
### Done reading OUTCAR file

############# calculate length quantum #############
length_quantum=[]
for i in range(3):
        tmp = 0
        for j in range(3):
                tmp = tmp + prim_mat[j][i]
        length_quantum.append(tmp)
length_quantum = np.array(length_quantum)
###############

##########  
print "   TOTAL      Ionic contribution \t:", ionic_contribution
print "   TOTAL electronic contribution \t:",ele_contribution
print "                   Voluem of cell\t:",omega
print "                Primitive Vector \t:",prim_mat
print "                  Length quantum \t:", length_quantum
total_pol = ionic_contribution + ele_contribution
print "             TOTAL Dipole Moment \t:", total_pol
quantum=[]
for i in range(3):
        quantum.append(round(total_pol[i]/length_quantum[i]))
quantum = np.array(quantum)
print "     Quantum X,Y and Z direction \t:",quantum
print "  Length Quantum along X,Y and Z \t:",quantum* length_quantum
polarization = -1*1602*( total_pol - (quantum* length_quantum) )/omega
print "    Components of   Polarization \t:",polarization
print "    Total Polarization  (uC/cm^2)\t:", np.linalg.norm(polarization)
