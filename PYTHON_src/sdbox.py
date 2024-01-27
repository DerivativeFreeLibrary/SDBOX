#============================================================================================
#    SDBOX - PYTHON implementation of a Derivative-Free algorithm for bound 
#    constrained optimization problems 
#    Copyright (C) 2017  G.Liuzzi
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#    S. Lucidi, M. Sciandrone. A Derivative-Free Algorithm for Bound Constrained Optimization, 
#    Computational Optimization and Applications, 21(2): 119-142 (2002)
#    DOI: 10.1023/A:1013735414984
#
#============================================================================================
import math
import problem
funct = problem.funct
setbounds = problem.setbounds
startp = problem.startp

#import problem1
#funct = problem1.funct
#setbounds = problem1.setbounds
#startp = problem1.startp

sqrt = math.sqrt

#----------------------------------------------------------
# function arresto
#----------------------------------------------------------
def arresto(alpha_tilde,tol,kfstop,fstop,fstoptol):
	n = len(alpha_tilde)
	if max(alpha_tilde) <= tol:
		return True
	else:
		if kfstop == n+1:
			fm = 0.0
			for i in range(1,kfstop+1):
				fm += fstop[i-1]

			fm /= float(kfstop)
			sq = 0.0
			for i in range(1,kfstop+1):
				sq += (fstop[i-1] - fm)**2

			sq = sqrt(sq / float(kfstop))
			if sq <= fstoptol:
				return True
	
	return False
#----------------------------------------------------------

#----------------------------------------------------------
# function linesearch
#----------------------------------------------------------
def linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,nf):
	
	#initializations
	delta = 0.5
	f_ext = fz
	alpha_tilde = min(alpha_max,alpha/delta)
	z[i] = x[i]+alpha_tilde*d[i] 
	fz  = funct(z)
	nf += 1
	while (alpha_tilde < alpha_max ) and (fz<=f-gamma*alpha_tilde**2):
		f_ext = fz
		alpha = alpha_tilde
		alpha_tilde = min(alpha_max,alpha/delta)
		z[i] = x[i]+alpha_tilde*d[i]
		fz  = funct(z)
		nf += 1
	return alpha, f_ext, nf
#----------------------------------------------------------

maxiter = 20000

lb, ub = setbounds()
x = startp()

n  = len(x)
d = [1.0]*n
fstop = [None]*(n+1)
alpha_tilde = [1.0]*n
nf = 0
gamma = 1.e-6	
theta = 0.5
tol = 1.e-6
fstoptol = 1.e-13

f = funct(x)
nf += 1
f_ext = f
ifstop = 1
kfstop = 1
fstop[ifstop-1] = f

FORMAT = "obj.f:%11.4e max alpha:%10.4e k:%6d nf:%6d i_coor:%3d"

# MAIN LOOP
for k in range(1,maxiter+1):

	z = [a for a in x]
	if arresto(alpha_tilde,tol,kfstop,fstop,fstoptol):
		break
	for i in range(n):
		print(FORMAT % (f, max(alpha_tilde),k,nf,i))
		flag_dir = True
		#--------------
		#STEP 1--> di
		#--------------
		f = f_ext
		if d[i] > 0.0:
			alpha_max = ub[i]-x[i]
		else:
			alpha_max = x[i]-lb[i]

		alpha = min(alpha_tilde[i],alpha_max)
		z[i] = x[i]+alpha*d[i]
		fz  = funct(z)
		nf += 1
		if (alpha > 0.0) and (fz <= f-gamma*alpha**2):
			#--------------
			#STEP 3(1)
			#--------------
			alpha, f_ext, nf = linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,nf)
			alpha_tilde[i] = alpha
			flag_dir = False

		#--------------
		#STEP 2--> -di
		#--------------
		if flag_dir:
			if d[i] > 0.0:
				alpha_max = x[i]-lb[i]
			else:
				alpha_max = ub[i]-x[i]

			alpha = min(alpha_tilde[i],alpha_max)
			z[i] = x[i]-alpha*d[i]
			fz  = funct(z)
			nf += 1
			if (alpha > 0.0) and (fz <= f-gamma*alpha**2):
				#--------------
				#STEP 3(2)
				#--------------
				d[i] = -d[i]
				alpha, f_ext, nf = linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,nf)
				alpha_tilde[i] = alpha
			else:
				alpha = 0.0
				alpha_tilde[i] = theta*min(alpha_tilde[i],alpha_max)

		#--------------
		#STEP 4
		#--------------
		if alpha > 0.0:
			ifstop = ifstop + 1
			if ifstop <= n+1:
				if kfstop < ifstop: 
					kfstop = ifstop
			else:
				ifstop = 1

			fstop[ifstop-1] = f_ext

		x[i] = x[i] + alpha*d[i]
		z[i] = x[i]
		f = f_ext

print(FORMAT % (f, max(alpha_tilde),k,nf,0))
