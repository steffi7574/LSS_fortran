      INTEGER n
      PARAMETER (n=1000)
      REAL*8 mu, dt,dtref

      REAL*8 x,v
      COMMON /primvar/ x(n),v(n)


      COMMON /const/ mu, dt,dtref
