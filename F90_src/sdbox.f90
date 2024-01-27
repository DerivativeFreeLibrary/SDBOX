!============================================================================================
!    SDBOX - FORTRAN90 implementation of a Derivative-Free algorithm for bound 
!    constrained optimization problems 
!    Copyright (C) 2011  G.Liuzzi, A. Risi
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!    S. Lucidi, M. Sciandrone. A Derivative-Free Algorithm for Bound Constrained Optimization, 
!    Computational Optimization and Applications, 21(2): 119-142 (2002)
!    DOI: 10.1023/A:1013735414984
!
!============================================================================================

PROGRAM sdbox
	IMPLICIT NONE
	INTEGER:: n, i, k, h
	INTEGER:: maxiter, nf
	LOGICAL:: flag_dir

	PARAMETER(maxiter=20000)

	DOUBLE PRECISION, ALLOCATABLE:: x(:), LB(:), UB(:), d(:), z(:), fstop(:), alpha_tilde(:)
	DOUBLE PRECISION:: f, fz, f_ext, alpha_max, alpha, tol
	DOUBLE PRECISION:: fstoptol, gamma, theta
	LOGICAL:: arresto
	INTEGER:: ifstop, kfstop

	! TIMING
	DOUBLE PRECISION:: tstart, tend
	
	call setdim(n)

	allocate(x(n), LB(n), UB(n), d(n), z(n), fstop(n+1), alpha_tilde(n))

	call startp(n,x)
	call setbounds(n,lb,ub)

	!INIZIALIZZAZIONI
	nf = 0
	alpha_tilde = 1.d0
	d = 1.d0
	gamma = 1.d-6	
	theta = 0.5d0
	tol = 1.d-6
	fstoptol = 1.d-6

	call funct(n,x,f)
	nf = nf+1
	f_ext = f
	ifstop = 1
	kfstop = 1
	fstop(ifstop) = f

	! CICLO PRINCIPALE
	do k=1, maxiter
		!WRITE(*,10) f, maxval(alpha_tilde),k,nf
		z = x
		if(arresto(n,alpha_tilde,tol,kfstop,fstop,fstoptol)) exit
		do i=1,n
			WRITE(*,10) f, maxval(alpha_tilde),k,nf,i
			flag_dir = .true.
			!--------------
			!STEP 1--> di
			!--------------
			f = f_ext
			if(d(i) > 0.d0) then
				alpha_max = UB(i)-x(i)
			else
				alpha_max = x(i)-LB(i)
			end if
			alpha = min(alpha_tilde(i),alpha_max)
			z(i) = x(i)+alpha*d(i)
			call funct(n,z,fz)
			nf = nf+1
			if((alpha > 0.d0).and.(fz <= f-gamma*alpha**2.d0)) then
				!--------------
				!STEP 3(1)
				!--------------
				call linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,f_ext,nf)
				alpha_tilde(i) = alpha
				flag_dir = .false.
			end if
			!--------------
			!STEP 2--> -di
			!--------------
			if(flag_dir) then
				if(d(i) > 0.d0) then
					alpha_max = x(i)-LB(i)
				else
					alpha_max = UB(i)-x(i)
				end if
				alpha = min(alpha_tilde(i),alpha_max)
				z(i) = x(i)-alpha*d(i)
				call funct(n,z,fz)
				nf = nf+1
				if((alpha > 0.d0).and.(fz <= f-gamma*alpha**2.d0)) then
					!--------------
					!STEP 3(2)
					!--------------
					d(i) = -d(i)
					call linesearch(i,n,x,z,d,alpha_max,alpha,	& 
							gamma,f,fz,f_ext,nf)
					alpha_tilde(i) = alpha
				else
					alpha = 0.d0
					alpha_tilde(i) = theta*min(alpha_tilde(i),alpha_max)
				end if
			end if
			!--------------
			!STEP 4
			!--------------
			if(alpha > 0.d0) then
				ifstop = ifstop + 1
				if(ifstop <= n+1) then
					if(kfstop < ifstop) kfstop = ifstop
				else
					ifstop = 1
				endif
				fstop(ifstop) = f_ext
			endif
			x(i) = x(i) + alpha*d(i)
			z(i) = x(i)
			f = f_ext
			!write(*,*) 'x:',x
		end do
	end do   
	
	WRITE(*,10) f, maxval(alpha_tilde),k,nf

	10&
	FORMAT('obj.f:'es11.4,' max alpha:',d10.4,' k:',i6,' nf:',i6,' i_coor:',i3)
		
END PROGRAM sdbox

LOGICAL FUNCTION arresto(n,alpha_tilde,tol,kfstop,fstop,fstoptol)
	IMPLICIT NONE
	INTEGER:: n, kfstop, i
	DOUBLE PRECISION:: alpha_tilde(n), tol,	fstop(n+1), fstoptol
	DOUBLE PRECISION:: fm,sq
	IF(maxval(alpha_tilde) <= tol) THEN
		arresto = .true.
		RETURN
	ELSE
		IF(kfstop == n+1) THEN
			fm = 0.d0
			DO i = 1,kfstop
				fm = fm + fstop(i)
			ENDDO
			fm = fm / dble(kfstop)
			sq = 0.d0
			DO i = 1,kfstop
				sq = sq + (fstop(i) - fm)**2.d0
			ENDDO
			sq = DSQRT(sq / dble(kfstop))
			IF(sq <= fstoptol) THEN
				arresto = .true.
				RETURN
			ENDIF
		ENDIF
	ENDIF
	
	arresto = .false.
	RETURN

END FUNCTION arresto

SUBROUTINE linesearch(i,n,x,z,d,alpha_max,alpha,gamma,f,fz,f_ext,nf)
	IMPLICIT NONE
	INTEGER:: n, i, nf
	DOUBLE PRECISION:: x(n), z(n), f, fz, gamma, d(n),f_ext
	DOUBLE PRECISION:: alpha_max, alpha, delta, alpha_tilde
	
	!INIZIALIZZAZIONI
	delta = 0.5d0
	f_ext = fz
	alpha_tilde = min(alpha_max,alpha/delta)
	z(i) = x(i)+alpha_tilde*d(i) 
	call funct(n,z,fz)
	nf = nf+1
	do while((alpha_tilde < alpha_max ).and.(fz<=f-gamma*alpha_tilde**2.d0))
		f_ext = fz
		alpha = alpha_tilde
		alpha_tilde = min(alpha_max,alpha/delta)
		z(i) = x(i)+alpha_tilde*d(i)
		call funct(n,z,fz)
		nf = nf+1
	end do
END SUBROUTINE linesearch

