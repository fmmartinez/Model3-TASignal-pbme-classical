#!/usr/bin/python
import os
from sys import argv

if len(argv) == 1:
	print "you must provide name of simulation (folder name) as an argument"
	quit()

dataname = argv[1]

#30 directories plus the zero one
np = 31
#truncation for running a short part of the trajectory (top should be less than 31)
top = 31

d = []

#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = dataname + '/e' + str(j) + '-sec0' + str(i)
		d.append(dir)

	for i in range(10,np):
		dir = dataname + '/e' + str(j) + '-sec' + str(i)
		d.append(dir)

currentpath = os.getcwd()
for j in range(0,8):
	for i in range(0,top):
		zz = i + j*np 
		workingpath = currentpath + '/' + d[zz]
		os.chdir(workingpath)
		os.system('qsub submit.pbs')
