C################################################################
      PROGRAM vanderpol
C################################################################
      IMPLICIT NONE
      INTEGER I,J,K,L
      REAL*8 residuum
      include "var.inc"

      INTEGER fileid
      CHARACTER*160 filename



c initialize
      dtref=0.05
      dt=dtref
      mu=1.0
      
      DO I=1,n
         x(I)=1.0
         v(I)=1.0
      ENDDO

c compute first residuum
      CALL resid(residuum)
      write(*,*) 'Residuum 0', residuum

c first k-iteration
      CALL primal(-1)
      DO I=1,n
      ENDDO
      CALL resid(residuum)
      write(*,*) 'Residuum 1', residuum

c open residuum output file
      open(unit=11,file='resid',status='replace',position='append')


c start k-loop
      DO K=1,150

c       ... open files for output
        fileid=k+100
        write(filename,"(I3,A1)") fileid,'x'
        IF (mod(k,10).EQ.0) open(unit=fileid,file=filename,
     &            status='replace', position='append')

c        ... update primals
        CALL primal(k)
        CALL resid(residuum)
        write(*,*) K,residuum
        write(11,*) K,residuum


c       ... close files
        IF (mod(k,10).EQ.0) close(unit=fileid)
      ENDDO
c end of k-loop


      close(11)

      END



C###############################################################
      SUBROUTINE resid(residuum)
C###############################################################
      IMPLICIT NONE
      INTEGER I
      REAL*8 rxnorm,rvnorm
      REAL*8 residuum
      include "var.inc"
      REAL*8 rx(n),rv(n)
      
      DO I=2,n
        rx(1)=0.0
        rv(1)=0.0
        rx(I)= 1.0/dt*(x(I)-x(I-1)) - v(I)
        rv(I)=1.0/dt*(v(I)-v(I-1)) + x(I) - mu*(1.0-x(I)**2)*v(I)
      ENDDO
      CALL normsquare(rx,n,rxnorm)
      CALL normsquare(rv,n,rvnorm)
      residuum = sqrt(rxnorm+rvnorm)
      
      RETURN
      END

C###############################################################
      SUBROUTINE primal(k)

c     Integer k = k-loop iteration number
C###############################################################
      IMPLICIT NONE
      INTEGER I,L
      INTEGER k, fileid
      REAL*8 f1, f2
      REAL*8 Jac(2,2)
      include "var.inc"

      DO I=2,n
        IF(k.EQ.1) THEN   !initialize x,v for Newtonsolver
          x(I)=x(I-1)
          v(I)=v(I-1)
        ENDIF
c       ...RHS evaluation
        f1 = 1.0/dt*(x(I)-x(I-1)) - v(I)
        f2 = 1.0/dt*(v(I)-v(I-1)) + x(I) - mu*(1.0-x(I)**2)*v(I)
c        ...setting up Jacobian of RHS
        Jac(1,1) = 1.0/dt
        Jac(1,2) = -1.0
        Jac(2,1) = 1.0+2.0*mu*x(I)*v(I)
        Jac(2,2) = 1.0/dt-mu*(1.0-x(I)**2)
        CALL invert2by2(Jac)
c       ...update x(I) and v(I) with relaxed Newton step
        x(I) = x(I) - 0.95*( Jac(1,1)*f1 + Jac(1,2)*f2 )
        v(I) = v(I) - 0.95*( Jac(2,1)*f1 + Jac(2,2)*f2 )
       
        fileid=k+100
        IF (k.GT.1.AND.mod(k,10).EQ.0.0) write(fileid,*) I, x(I)
      ENDDO

      END


C###############################################################
      SUBROUTINE normsquare(a,M,nor)
C###############################################################
      IMPLICIT NONE
      INTEGER I,M
      REAL*8 a(M),nor
      
      nor=0.0
      DO I=2,M
        nor=nor+a(I)**2
      ENDDO
      
      RETURN
      END

C###############################################################
      SUBROUTINE invert2by2(A)
C###############################################################
      IMPLICIT NONE
      REAL*8 A(2,2)
      REAL*8 det,temp
      INTEGER L

      det = A(1,1)*A(2,2)-A(1,2)*A(2,1)
      IF ((det.LE.1e-14).AND.(det.GE.-1e-14)) THEN
          write(*,*) "ERROR ! Jacobian not invertible"
          stop
      ELSE
        det=1.0/det
        temp = A(1,1)
        A(1,1) = det*A(2,2)
        A(2,2) = det*temp
        A(1,2) = -det*A(1,2)
        A(2,1) = -det*A(2,1)
      ENDIF

      END

