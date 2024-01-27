subroutine setdim(n)
	implicit none
	integer		:: n

	n = 4
	return
end subroutine setdim

SUBROUTINE setbounds(n,A,B)
	IMPLICIT none
	integer		:: n 
	real*8		:: A(n), B(n)	
	integer		:: i

	DO I=1,n
		A(I)=0.d0
		B(I)=10.d0
	ENDDO
	RETURN
END SUBROUTINE setbounds

SUBROUTINE startp(n,X)
	IMPLICIT none
	integer		:: n 
	real*8		:: X(n)	
	integer		:: i

	DO I=1,n
		X(I)=5.d0
	ENDDO
	RETURN
END SUBROUTINE startp

SUBROUTINE FUNCT(N,X,F)

      IMPLICIT NONE

      INTEGER          :: N
      DOUBLE PRECISION :: X(N), F

      DOUBLE PRECISION :: A(10,4), C(10), FA
      INTEGER          :: I, J

      DO I=1,4
         A(1,I)=4.D0
         A(2,I)=1.D0
         A(3,I)=8.D0
         A(4,I)=6.D0
      END DO
      DO I=1,2
         A(5,2*(I-1)+1)=3.D0
         A(5,2*I)=7.D0
         A(6,2*(I-1)+1)=2.D0
         A(6,2*I)=9.D0
         A(7,I)=5.D0
         A(7,I+2)=3.D0
         A(8,2*(I-1)+1)=8.D0
         A(8,2*I)=1.D0
         A(9,2*(I-1)+1)=6.D0
         A(9,2*I)=2.D0
         A(10,2*(I-1)+1)=7.D0
         A(10,2*I)=3.6D0
      END DO

      C(1)=0.1D0
      C(2)=0.2D0
      C(3)=0.2D0
      C(4)=0.4D0
      C(5)=0.4D0
      C(6)=0.6D0
      C(7)=0.3D0
      C(8)=0.7D0
      C(9)=0.5D0
      C(10)=0.5D0

      F  = 0.0D0
      FA = 0.0D0

      DO I=1,10
         DO J=1,4
            FA = FA +(X(J)-A(I,J))**2
         END DO
         IF ((FA+C(I)).EQ.0.D0) THEN
            F=1.D25
            RETURN
         ENDIF
         F  = F -1.0D0/(FA+C(I))
         FA = 0.0D0
      END DO

      RETURN

END SUBROUTINE FUNCT

