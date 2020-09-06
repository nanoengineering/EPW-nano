#!/bin/sh
####################################################################
# Copyright Eyvaz Isaev
#
#
# Department of Physics, Chemistry, and Biophysics (IFM), 
# Linkoping University, Sweden  
#
# Theoretical Physics Department,
# Moscow State Institute of Steel and Alloys, Russia
# (Technological University)
# 
# isaev@ifm.liu.se, eyvaz_isaev@yahoo.com 
#
# set environment variables
#
. ../environment_variables
#
# As input parameters you have to specify (there are no default parameters):
#
# Mandatory parameters:
#
# Atoms - specify atoms involved in your calculations, you can specify Fe1, Fe1, Cr1,Cr2, etc.
#         N.B.! in the same order as you did it in your self-consistent input file.
# delta_e - frequency step (in cm{-1})
#
# "Temperature" file which contains
# T_start, T_end, T_delta
# initial T, final T, and delta T for calculation of the phonon contribution (F_vib) to 
# the free energy, the heat capacity (C_v), phonon entropy (S_vib), and phonon internal 
# energy (E_int)  
#
# T_list - list of temperatures for which the Debye temperature has to be calculated.
#          For lower T (do not specify T<10K !!!) calculations require more computer time,
#          starting from 10K is much faster 
#
# Optional parameters:
# Sysname - specify your system name you have considered, 
#           required to see your system name in output files
# SysInfo - Any information on your system (volume, pressure, etc.) will be reflected in the head part 
#           of output file 
#
# Input parameters for Quasiharonic calculations (fqha.in)
# PHDOS.out - Total phonon DOS, you might specify any partial phonon DOS name but only 
#             total phonon DOS is used for this purpose
# Sysname.QHA.out - Name of output file for QHA calculations
# T_start, T_end, T_step - starting T, the highest T used in QHA calculations, and T step
#
#
# To start Phonon DOS calculatons you have to edit matdin.init file (See below) to specify 
#          force constants matrix, and atomic masses
# You need also "ttrinp" file to show Brillouin Zone:
#          See ttrinp file for explanation how  you can manage it
# There are some ttrinp files for a number of popular (FCC, BCC, Simple Cubic, and HCP) latticies
# For these crystal structures you need no additional information and ttrinp file added 
# automatically.
##############################################################################  

##############################################################################  
# Optional parameters, any information specific for the system studied 
# 
SysInfo='InSb'

# Mandatory parameters
# Specify SystemName and Force Constants matrix

Sysname='InSb'
FC_file='InSb666.fc'

#
# Specify lattice type (used to create ttrinp file). It should be the same as in scf.in file 
# Specify atoms in the unit cell as they specified in scf.in file
# Specify atomic masses for these atoms in the same order as in scf.in
# Specify the frequency step (delta_e) as well, but 0.75 is a good choice

ibrav=2
atoms="In   Sb  "
mass="114.818  121.76 "
delta_e=0.75

# Edit ONLY amass parameters
# Please do not change flfrq='frequency' line
# leave asr (acoustic sum rule) and flfrc lines

cat >matdyn.init <<EOF
&input
    amass(1)=114.818,
    amass(2)=121.76,
    asr='crystal',
    flfrc='$FC_file',
    flfrq='frequency'
/
EOF

#
# In most cases there is no need to edit files listed below, but if you like ...
#

# Temperature range for thermodynamic properties
# T_start, T_end, T_step for QHA calculations

cat > Temperature <<EOF
5 500 5
EOF

# Debye Temperature calculations
# Phonon DOS filename (total phonon DOS, not projected), leave it as PHDOS.out
# accuracy (limited 1.d-5, more accuracy is not required )
# Low_Temp_start, Low_Temp_end, and Low_Temp_step for Low Temperature  limit, up to 15-30K
# Hihg temperature and T_step for HT limit

cat >T_Debye.in <<EOF
PHDOS.out
0.0001
 3 15 3
 500 10
EOF



#############################################################################
#
# Copy appropriate tetrahedra file to ttrinp 
#

case $ibrav  in

0)   echo "You should know about the symmetry of a system you study"
     exit ;;
1)   cp $QHA_DIR/tetrahedra/ttrinp_sc ./ttrinp  ;;
2)   cp $QHA_DIR/tetrahedra/ttrinp_fcc ./ttrinp ;;
3)   cp $QHA_DIR/tetrahedra/ttrinp_bcc ./ttrinp ;;
4)   c2a1=`head -1 $FC_file | cut -c 36-44`
     c2a=`echo "scale=8;(1/$c2a1)/2" | bc -l`
     sed 's/X/'$c2a'/g' $QHA_DIR/tetrahedra/ttrinp_hcp > ./ttrinp ;;
5)   echo "Trigonal R: not implemented yet"
     exit ;;
6)   c2a1=`head -1 $FC_file | cut -c 36-44`
     c2a=`echo "scale=8;(1/$c2a1)/2" | bc -l`
     echo $c2a
     sed 's/X/'$c2a'/g' $QHA_DIR/tetrahedra/ttrinp_stetra > ./ttrinp ;; 
7)   echo "Running script instruction:"
     echo "See instructions in $QHA_DIR/tetrahedra/ttrinp_bct file to setup tetrahedra vertecies for c/a<1"
     echo "Then comment exit line by #"
     echo "And uncomment the next line"
     exit;;
#    cp $QHA_DIR/tetrahedra/ttrinp_bct ./ttrinp ;;
8)   b2a1=`head -1 $FC_file | cut -c 25-33`
     c2a1=`head -1 $FC_file | cut -c 36-44`
     b2a=`echo "scale=8;(1/$b2a1)/2" |bc -l`
     c2a=`echo "scale=8;(1/$c2a1)/2" |bc -l`
     sed 's/YY/'$b2a'/g' $QHA_DIR/tetrahedra/ttrinp_ortho_simple > ttrinp1
     sed 's/ZZ/'$c2a'/g' ttrinp1 > ./ttrinp
     rm -f ttrinp1;;
9)   echo "Orthorhombic base centered: not implemented yet"
     exit ;;
10)  echo "Orthorhombic face centered: not implemented yet"
     exit ;;
11)  echo "Orthorhombic body centered: not implemented yet"
     exit ;;
12)  echo "Monoclinic P: not implemented yet"
     exit ;;
13)  echo "Monoclinic base centered: not implemented yet"
     exit ;;
14)  echo "Triclinic: not implemented yet"
     exit ;;
esac 

############################################################################
# Below run commands 

# Generate q-points
$QHA_DIR/bin/tetra.x 

cp matdyn.init matdyn.init.tmp
cat  kpts_out >> matdyn.init.tmp
echo  EOF >>matdyn.init.tmp

mv  matdyn.init.tmp matdyn.in

echo ' Recalculating omega(q) from C(R)'
$BIN_DIR/matdyn.x < matdyn.in > matdyn.out

nmodes=`head -1 frequency | cut -c 13-16 `
nkpt=`head -1 frequency | cut -c 23-26 `
natoms=`echo "scale=0; $nmodes/3" | bc -l`
#
cat >phdos1.in <<EOF
$nkpt  $nmodes
$atoms
EOF

# Calculate partial phonon DOS
$QHA_DIR/bin/Partial_phonon_DOS.x < phdos1.in 

#rm -f phdos1.in

cat >phdos.in <<EOF
$delta_e
$atoms
EOF

# Calculate total  phonon DOS and atom projected phonon DOS
$QHA_DIR/bin/phonon_dos.x <frequency 

#rm -f phdos.in

# Remove NaN from  all phonon DOS files

cat >name <<EOF
PHDOS.out
EOF

cp PHDOS.out PHDOS.out.copy

$QHA_DIR/bin/Ghost_DOS.x <name  >out 

mv out PHDOS.out
#rm -f name

# Atomic related properties

cat >atom_info <<EOF
$natoms
$atoms
$mass
EOF

for Atom in $atoms

do 

cat > atom_name <<EOF
$Atom
EOF


$QHA_DIR/bin/atom_info.x < atom_info > atom_mass

cat >name <<EOF
projected_DOS.$Atom
EOF

$QHA_DIR/bin/Ghost_DOS.x <name  >out

mv out projected_DOS.$Atom

cp projected_DOS.$Atom projected.DOS

echo "# $Sysname $Atom  $SysInfo" >>Thermodynamics.$Atom

$QHA_DIR/bin/Atom_projected_properties.x >>Thermodynamics.$Atom

# 
# Mean Square Displacement calculations for each atoms 

cat name Temperature atom_mass > displacement.in

$QHA_DIR/bin/Mean_square_displacement.x < displacement.in

mv Displacements Displacements.$Atom

done

# Debye Temperature calculations

$QHA_DIR/bin/Debye.x >> Theta_D

#rm -f T_Debey.in

# Finally, thermodynamic properties
#
# Parameters required for QHA calculations
# Total Phonon DOS file
# output file for C_V, S, Internal energy
#
cat >fqha.in <<EOF
PHDOS.out
$Sysname.QHA.out
EOF

cat Temperature >> fqha.in

$QHA_DIR/bin/F_QHA.x <fqha.in 

#rm -f fqha.in

echo 'Phonon DOS and Quasiharmonic calculations have finished.' 
echo 'Now you can analyse these data using Gnuplot or xmgrace'
echo 'Enjoy!'


