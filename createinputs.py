#!/usr/bin/python
import os
import shutil

#run value
deltahere = 0.0

#new values
dt    = '5d-5'
nmds  = '954'
step  = '64'

nmcs = '100000'

#original values
#double check that dt, nmds and step are related with these!!!
#for example halving dt doubles nmds and step
#nmds0 = '477'
#dt0   = '1.d-4'
#step0 = '64'

#number of runs per signal, n+1 because 0 is included
np = 61

#if np doubles this number halves, as well as step above
step2 = 0.02

#number of basis functions per center
bg = 1
bb = 1
bd = 1

#pulses stuff
e0 = 0.9
e1 = 0.09

tau1 = '0.045'
omega1 = '260'
tau2 = '0.045'
omega2 = '260'

#walltime for simulations
wallt = '1:00:00'

d = []

#make root directory
rod = 'TA%s%s%s' % (str(bg),str(bb),str(bd))
#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = rod + '/e%s-sec0%s' % (str(j),str(i))
		d.append(dir)

	for i in range(10,np):
		dir = rod + '/e%s-sec%s' % (str(j),str(i))
		d.append(dir)

#Generate Folders if not present
if not(os.path.exists(rod)):
	os.mkdir(rod)

for i in range(0,np*8):
   if not(os.path.exists(d[i])):
      os.mkdir(d[i])

#############
#	Generation of submit.pbs, map.in and intesity.in
for k in range(0,8):
	if (k == 0):
		p1 = 0
		p2 = 0
		p3 = 0
		p0name = rod + '_000-'
	if (k == 1):
		p1 = 0
		p2 = 1
		p3 = 0
		p0name = rod + '_010-'
	if (k == 2):
		p1 = 0
		p2 = 0
		p3 = 1
		p0name = rod + '_001-'
	if (k == 3):
		p1 = 0
		p2 = 1
		p3 = 1
		p0name = rod + '_011-'
	if (k == 4):
		p1 = 1
		p2 = 0
		p3 = 0
		p0name = rod + '_100-'
	if (k == 5):
		p1 = 1
		p2 = 1
		p3 = 0
		p0name = rod + '_110-'
	if (k == 6):
		p1 = 1
		p2 = 0
		p3 = 1
		p0name = rod + '_101-'
	if (k == 7):
		p1 = 1
		p2 = 1
		p3 = 1
		p0name = rod + '_111-'
	
	for i in range(0,np):
		ii = i + k*np
		mapfile = open('./' + d[ii] + '/map.in','w')
		intfile = open('./' + d[ii] + '/intensity.in','w')
		subfile = open('./' + d[ii] + '/submit.pbs','w')
#		print 'cp map ../'+d[ii]+'/.'

#write each map.in file
		l='''Np	DELTA	KONDO	NOSC	OME_MAX
3	%s	0.09	20	50
NMCS	NMDS	seed	DT	LUMDA_D
%s	%s	90	%s	10
Eg	Eb	Ed	mu	E0	E1	beta	vib_omega
0	240	240	1	%s	%s	0.24	37.7
TAU1	OMEGA1	TAU2	OMEGA2	time3	step1	step2
%s	%s	%s	%s	0.15	%s	%s
BATH(0:B EQ 1:T EQ)	INIT	NFILE
0	3	%s
basisfg	b	d
%s	%s	%s
pi	pj	pk
%s	%s	%s''' % (str(deltahere),nmcs,nmds,dt,e0,e1,tau1,omega1,tau2,omega2,
					step,str(step2),str(i),str(bg),str(bb),str(bd),str(p1),str(p2),str(p3))

		mapfile.write(l)

		nmdsstep = str(int(nmds) + i*int(step))
		delay = str(0.15 + i*step2)

#write each intensity.in file
		m='''E1	NMDS	tau	omega	delay
%s	%s	%s	%s	%s''' % (str(e1),nmdsstep,tau2,omega2,delay)
		
		intfile.write(m)
	
#write each submit.pbs
		
		n = '''#!/bin/bash -l
#PBS -S /bin/bash
#PBS -r n
#PBS -N %s%s
#PBS -l walltime=%s
#PBS -l procs=1
#PBS -l mem=256mb

cd $PBS_O_WORKDIR
time ./a.out < map.in''' %(p0name,str(i),wallt)
		
		subfile.write(n)
		
		#
		intfile.close()
		mapfile.close()
		subfile.close()

#copy executables
shutil.copy2('int.out',rod)
for i in range(0,8*np):
	shutil.copy2('a.out',d[i])
