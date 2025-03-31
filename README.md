# silkmdsimulation

#####Download LAMMPS (GPU version of bump included)#####
git clone https://github.com/mppotter/lammps.git

#####create a silk chain unit#####
#set up simulation box

variable        dimx equal 16
variable        dimy equal 7
variable        dimz equal 7
region		box block  0  ${dimx}  2  ${dimy}  2  ${dimz}
create_box	9 box


#make custumized lattice for a hard segment
variable        rho2 equal 18./(2)^3 
lattice	custom ${rho2} & 
                basis 0.01 0.51 0.51 basis 0.51 0.51 0.51 &
                basis 0.01 0.26 0.51 &
                basis 0.01 0.51 0.26 &
                basis 0.01 0.76 0.51 &
                basis 0.01 0.51 0.76 &
                basis 0.51 0.26 0.51 &
                basis 0.51 0.51 0.26 &
                basis 0.51 0.76 0.51 &
                basis 0.51 0.51 0.76 &
                basis 0.01 0.26 0.26 &
                basis 0.01 0.76 0.26 &
                basis 0.01 0.76 0.76 &
                basis 0.01 0.26 0.76 &
                basis 0.51 0.26 0.26 &
                basis 0.51 0.76 0.26 &
                basis 0.51 0.76 0.76 &
                basis 0.51 0.26 0.76 &

#create a hard segment in the middle of a chain unit
region		middle block  6 10 5 6 5 6   #H8S24 25% hard segment

#H4S28 (12.5%) 7 9 5 6 5 6 
#H12S20 (37.5%) 5 11 5 6 5 6 
#H16S16 (50%) 4 12 5 6 5 6 

create_atoms	9 region middle basis 1 1 basis 2 2 &
		basis 3 3 &
		basis 4 4 &
		basis 5 3 &
		basis 6 4 &
		basis 7 3 &
		basis 8 4 &
		basis 9 3 &
		basis 10 4 &
		basis 11 6 &
		basis 12 6 &
		basis 13 6 &
		basis 14 6 &
		basis 15 6 &
		basis 16 6 &
		basis 17 6 &
		basis 18 6 &

#create a soft segment
region		middle5 block  0 6 5 6 5 6
create_atoms	9 region middle5 basis 1 7 basis 2 8

region		middle2 block  10 16 5 6 5 6
create_atoms	9 region middle2 basis 1 7 basis 2 8

#create solvent particles (50wt%)
region		middle5 block  0 6 5 6 5 6
create_atoms	9 region middle5 basis 1 7 basis 2 8

region		middle2 block  10 16 5 6 5 6
create_atoms	9 region middle2 basis 1 7 basis 2 8

#mass of particles (soft segment = 1.0 / hard segment = 1.0 (0.1*8outerparticles+0.2centerparticles) / solvent = 0.35)
mass		* 1.0
mass		1 0.2
mass		2 0.2
mass		3 0.1
mass		4 0.1
mass		5 0.35

mass		6 0.1
mass		7 1.0
mass    8 1.0

#group particles
group           patch4 type 4 
group		patch1 type 1
group		patch2 type 2
group		patch3 type 3
group		patch6 type 6
group           patch5 type 5
group           patch7 type 7
group           patch8 type 8

#temperature setup
variable	t equal 0.2
variable	tf equal 0.2
velocity	all create $t 55

######pair potentials#####
#for initial setup, no secondary interactions are turned on
pair_style	lj/cut/bump 2.8 #cutoff 2.8
pair_coeff      * * 0.02 $v 1.0 1 0 $u

variable	c equal 1.8
variable	s equal $c/2^(1./6.)
variable	u equal 1.0
variable	v equal $u/2^(1./6.)

variable	p equal 1.0
variable	q equal $p/2^(1./6.)

variable	h equal 1.0
variable	k equal $h/2^(1./6.)

variable	h2 equal 1.0
variable	k2 equal ${h2}/2^(1./6.)

variable	h3 equal sqrt(2.0) # 1.414
variable	k3 equal ${h3}/2^(1./6.)

variable        b equal 0.2
variable	b1 equal 0.14

variable        c3 equal 2.0
variable        s3 equal ${c3}/2^(1./6.)

variable        c5 equal 1.6
variable        s5 equal ${c5}/2^(1./6.)

variable        h4 equal 0.5
variable        k4 equal ${h4}/2^(1./6.)

variable        h5 equal 0.25
variable        k5 equal ${h5}/2^(1./6.)

variable        c4 equal 0.707 
variable        s4 equal ${c4}/2^(1./6.)

variable        c6 equal 2.5
variable        s6 equal ${c6}/2^(1./6.)

variable        c7 equal 2.8
variable        s7 equal ${c7}/2^(1./6.)

#backbone
#in hard segment
pair_coeff	1 1 5.0 $s  1.0 1 0 1.8 #1.8 only repulsion 
pair_coeff	2 2 5.0 $s  1.0 1 0 1.8 #1.8 only repulsion
pair_coeff	1 2 5.0 ${k2} 1.2 1 1.0 1.8 #1 strong bond  

#connect soft segments in outside and inside hard segments
pair_coeff      2 7 20.0 ${k2} 1.2 1 1.0 1.8 #1 strong bond 
pair_coeff      1 8 20.0 ${k2} 1.2 1 1.0 1.8 #1 strong bond 
pair_coeff      1 7 20.0 $s  1.0 1 0 1.8 #1.8 only repulsion
pair_coeff      2 8 20.0 $s  1.0 1 0 1.8 #1.8 only repulsion

#soft segment (distance between soft/soft = 2.2)
pair_coeff      7 7 20.0 $s  1.0 1 0 1.8 #1.8 only repulsion
pair_coeff      8 8 20.0 $s  1.0 1 0 1.8 #1.8 only repulsion
pair_coeff      7 8 20.0 ${k2} 1.2 1 1.0 1.8 #1 strong bond 

#intramolecular interactions in hard segments
pair_coeff	1 3 20.0 ${k4} 1.2 1 2.0 1.118 #0.5 strong bond
pair_coeff	2 3 20.0 ${k4} 1.2 1 2.0 1.118 #0.5 strong bond
pair_coeff	1 4 20.0 ${k4} 1.2 1 2.0 1.118 #0.5 strong bond
pair_coeff	2 4 20.0 ${k4} 1.2 1 2.0 1.118 #0.5 strong bond
pair_coeff      1*2 6 20.0 ${s4} 1.2 1 2.0 1.225 #0.707 strong bond

pair_coeff      3 4 20.0 ${s4} 1.2 1 2.0 1.225 #0.707 weak bond
pair_coeff      3*4 6 20.0 ${k4} 1.2 1 2.0 1.118 #0.5 strong bond
#pair_coeff      3*4 6 0.01 ${k5} 1.2 1 0.0 0.25 #0.25 only repulsion

pair_coeff	3 3 5.0 ${k2} 1.2 1 1.0 1.414 #1 strong bond
pair_coeff	4 4 5.0 ${k2} 1.2 1 1.0 1.414 #1 strong bond

pair_coeff      6 6 5.0 ${k2} 1.2 1 1.0 1.414 #1 strong bond

#intermolecular interactions between soft and hard segments

pair_coeff      3*4 7 0.01 ${k4} 1.2 1 0.0 0.5 #0.5 only repulsion
pair_coeff      3*4 8 0.01 ${k4} 1.2 1 0.0 0.5 #0.5 only repulsion
pair_coeff      6 7 0.01 ${k4} 1.2 1 0.0 0.5 #0.5 only repulsion
pair_coeff      6 8 0.01 ${k4} 1.2 1 0.0 0.5 #0.5 only repulsion

#interaction with solvent particles
variable        c8 equal 2.0
variable        s8 equal ${c8}/2^(1./6.)

variable        c9 equal 2.2
variable        s9 equal ${c9}/2^(1./6.)

variable        c10 equal 0.3
variable        s10 equal ${c10}/2^(1./6.)

pair_coeff      1*2 5 0.05 ${s9} 1.2 1 0.0 2.2 # 2.2 only repulsion

pair_coeff      3*4 5 0.01 ${s10} 1.2 1 0.0 0.3 #only repulsion
pair_coeff      5 6 0.01 ${s10} 1.2 1 0.0 0.3 #only repulsion

pair_coeff      5 5 0.05 ${s9} 1.2 1 0.0 2.2 # 2.2 only repulsion

pair_coeff      5 7 0.05 ${s8} 1.2 1 0.0 2.0 # 2.0 only repulsion
pair_coeff      5 8 0.05 ${s8} 1.2 1 0.0 2.0 # 2.0 only repulsion


pair_modify	tail no
pair_modify	shift yes
neighbor	0.5 bin
neigh_modify	every 1 delay 0 check yes 

timestep	0.001
thermo		100
thermo_style    custom step temp pe ke etotal pxx pyy pzz lx ly lz press vol  

#####create a silk chain#####
variable        nrepeat equal 16 #chain length
replicate       ${nrepeat} 1 1
region  end block  0 2 INF INF INF INF units box
delete_atoms    region end

##########################step1: Initial setup including relaxation#################################
#increase the box
variable        newsized equal ly*16
fix             2 all deform 1 y final 0 ${newsized} z final 0 ${newsized} units box
run             100000

#deform to x axis
unfix           2
variable        lengthx equal 50*${nrepeat}/4.0
fix             2 all deform 1 x final 0 ${lengthx} units box
run             100000
unfix           2

#make the simulation box cubic
variable        newsize equal (lx*ly*lz)^(1./3.)
fix             2 all deform 1 x final 0 ${newsize} y final 0 ${newsize} z final 0 ${newsize} units box
run             100000
unfix           2

#decide # of chain in a simulation box
replicate       72 4 4   #replicate into 1152 chains
velocity        all create ${tf} 9936    #add velocity

#high temp relaxation
unfix           1
fix             1 all nvt temp ${tf} ${tf} 0.1
run             700000

#Increasing density 1000 per 1unit
variable        newsizerx equal lx/(320^(1./3.)) 
variable        newsizery equal ly/(320^(1./3.))
variable        newsizerz equal lz/(320^(1./3.))

fix             2 all deform 1 x final 0 ${newsizerx} y final 0 ${newsizery} z final 0 ${newsizerz} units box
run             500000
unfix           2

#write_data      s0.data

#cooling
variable        ct equal 0.04
unfix           1
fix             1 all nvt temp ${tf} ${ct} 0.1 #cooling
run             200000
#relaxation by maintaining temp
unfix           1
fix             1 all nvt temp ${ct} ${ct} 0.01
run             300000
unfix           1

write_data      s1.data 
quit

#########################step2: assembly simulation#############################
#turn on hydrophobic interaction of solvent particles
variable        c8 equal 2.0
variable        s8 equal ${c8}/2^(1./6.)

variable        c9 equal 2.2
variable        s9 equal ${c9}/2^(1./6.)

variable        c10 equal 0.3
variable        s10 equal ${c10}/2^(1./6.)

pair_coeff      1*2 5 0.05 ${s9} 1.2 1 0.0 2.2 # 2.2 only repulsion

pair_coeff      3*4 5 0.01 ${s10} 1.2 1 0.0 0.3 #only repulsion
pair_coeff      5 6 0.01 ${s10} 1.2 1 0.0 0.3 #only repulsion

#turn on 33bond and 44bond
pair_coeff      3 3 2.5 ${k2} 1.1 0.415 2.0 2.8 #1 strong bond
pair_coeff      4 4 2.5 ${k2} 1.1 0.475 2.0 2.8 #1 strong bond

#solvent behavior
pair_coeff      5 5 0.05 ${s9} 1.2 1 0.0 2.2 # 2.2 only repulsion

pair_coeff      5 7 0.05 ${s8} 1.2 1 0.0 2.8 # 2.0 attraction
pair_coeff      5 8 0.05 ${s8} 1.2 1 0.0 2.8 # 2.0 attraction


#relaxation 
fix		1 all nvt temp ${tf} ${tf} 0.01
run		500000 #pre-assemble
unfix           1
#P change into 0.001
variable        pressurenow equal press
fix             1 all npt temp ${tf} ${tf} 0.01 iso ${pressurenow} 0.001 10 
run             100000
unfix           1

#the first compression at low P (P = 0.001)
fix             1 all npt temp ${tf} ${tf} 0.01 iso 0.001 0.001 10

#dehydration (from 50wt% to 5wt%)
variable a loop 400
label loop
                group g5 type 5 5

variable        n5 equal count(g5)
variable        rmn equal 3925 #determine the final water concentration, need to calculate
variable        fraction equal ${rmn}/${n5}
                set group g5 type/fraction 9 ${fraction} 1168

                group g9 type 9 9
                delete_atoms group g9
                group g5 delete
                group g9 delete
run             5000
next a
jump s_hydro2 loop
unfix           1

############turn on Hbond in soft segments###################
pair_coeff      7 8 10.0 ${k2} 1.2 0.79 1.0 2.8 #1 strong bond/2.2 attraction

#P change into 1
variable        pressurenow equal press
fix             1 all npt temp ${tf} ${tf} 0.01 iso ${pressurenow} 1 10
run             100000
unfix           1

#the second compression at high P (P = 1)
fix		1 all npt temp ${tf} ${tf} 0.01 iso 1 1 10
run		300000
unfix           1

#decompression
variable       finalpress equal press
fix             1 all npt temp ${tf} ${tf} 0.01 iso ${finalpress} 0.0 10
run            200000
unfix          1

#maintain pressure_relaxation
fix             1 all npt temp ${tf} ${tf} 0.01 x 0.0 0.0 5.0 y 0.0 0.0 5.0 z 0.0 0.0 5.0
run            200000
unfix          1

write_data    resultingsilk.data


################step3: uni-axial tensile test##################################
