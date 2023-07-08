# naively written mc code for lj particles
# not cross-checked
# python python_MC.py 'nsteps' 'nparticles'
# e.g. python python_MC.py 100 125

import random
import numpy as np
import math
import subroutine
import sys
#from inputs import * 

print ( "total arguments passed", len(sys.argv))
print ( 'script name=',sys.argv[0])
print ('steps = ' ,sys.argv[1] )
print ('particles = ', sys.argv[2] )

#https://www.geeksforgeeks.org/command-line-arguments-in-python/



#***************************input and conversion
sigma=1.0
epsilon=1.0
#nstep=2000
nstep = int (sys.argv[1])
#n=125
n = int (sys.argv[2])
temp=1.2
beta=1/temp
dens=0.7             
box=(n/dens)**(1/3)
boxinv=1/box
rcut=(box/2)-.3      
rmin=0.7
drmax=.45
nadjst=10
histstep=0.0

max_bin=1000                         
delr=box/(2*max_bin)

file2 = open("potential.txt","w")
ff = open("trajectory.xyz","w")

#**************	simple cubic configuration => will eventually melt depending on t and rho*****************************

rx,ry,rz = np.zeros (n),np.zeros (n),np.zeros (n)
rxnew,rynew,rznew = np.zeros (n),np.zeros (n),np.zeros (n)

ll = n**(1/3)
print ( 'cubic box size = ',ll, 'simulation box size = ', box )

if ll %2 ==0:
	l1 = int (ll/2.0)
	l2 = l1
if ll %2 != 0:
	l1 = int (ll/2.0)
	l2 = l1+1 


jj = 0
for i in range (-l1,l2):
	for j in range (-l1,l2):
		for k in range (-l1,l2):
			rx[jj],ry[jj],rz[jj] = i, j, k
			jj = jj +1

if (jj != n):
	print ('something is wrong', jj, n)
	print ('put particle number cube of some odd numbers like 27,125....')

#************ core monte carlo part *******************

for step in range(0,nstep):
	naccept = 0.0
	ntrial = 0.0
	
	for i in range (0,n):
		v = 0.0
		en_return=subroutine.energy(n,rx[i],ry[i],rz[i],i,rcut,sigma,box,rx,ry,rz)
		vold=en_return
		
		
		rxnew[i] = rx[i] + ( 2.0 * random.random()- 1.0 ) * drmax
		rynew[i] = ry[i] + ( 2.0 * random.random()- 1.0 ) * drmax
		rznew[i] = rz[i] + ( 2.0 * random.random()- 1.0 ) * drmax

		rxnew[i] = rxnew[i] - round ( rxnew[i]*boxinv )*box		
		rynew[i] = rynew[i] - round ( rynew[i]*boxinv )*box
		rznew[i] = rznew[i] - round ( rznew[i]*boxinv )*box
		
		ntrial = ntrial + 1.0
		
		#[overlap checking]
		noverlap=subroutine.check(i,n,box,boxinv,rxnew[i],rynew[i],rznew[i],rmin,rx,ry,rz)
		
		if(noverlap==(n-1)):
			en_return=subroutine.energy(n,rxnew[i],rynew[i],rznew[i],i,rcut,sigma,box,rx,ry,rz)
			vnew=en_return
			
			deltv = (vnew-vold)
			deltvb = beta*(vnew-vold)
			if (deltvb < 75) :
				if(deltv < 0.0):
					v , naccept = v+deltv, naccept + 1.0
					rx[i],ry[i], rz[i] = rxnew[i], rynew[i], rznew[i]
					
				
				elif(math.exp(-deltvb) > random.random()):
					v , naccept = v+deltv, naccept + 1.0
					rx[i],ry[i], rz[i] = rxnew[i], rynew[i], rznew[i]
					
	if (step%nadjst==0):
		if ((naccept/ntrial) > .50):
			drmax=1.05*drmax
		else:
			drmax=.95*drmax

	if (step%5==0):
		vtotal=subroutine.sumup(rx,ry,rz,box,n,rcut,rmin,sigma)
		print(step,"%2.3f"%drmax,naccept/ntrial,"%0.4f"%vtotal,v)
		print (step, (vtotal/n), (v/n),sep='\t\t',file=file2)
		
	#*********** storing configuration *****************
	
	if (step%5==0):
		print (n,file=ff)
		print (n,step, file=ff)
		for i in range (n):
			rxstore = rx[i]-box*round(rx[i]*boxinv)
			rystore = ry[i]-box*round(ry[i]*boxinv)
			rzstore = rz[i]-box*round(rz[i]*boxinv)
			print ("c", rxstore, rystore, rzstore,file=ff)

file2.close()
ff.close()	
		
