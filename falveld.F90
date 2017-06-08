REAL(KIND=8) FUNCTION falveld(d,nu,s,g,p,csf)
  IMPLICIT NONE
! Computes particle fall velocity according to W. Dietrich, WRR 18(6),
! 1615-1626, 1982. falveld in m/s.
! Input parameters:
!   d    particle diameter (m)
!   nu   kinematic viscosity (m^2/s)
!   s    specific gravity
!   g    acceleration due to gravity (m/s^2)
!   p    Powers roundness factor
!   csf  Corey shape factor (must be >= 0.2)
!
! NOTE: for a sphere, p=6.0, csf=1.0; for natural sediments csf=0.7.
! Function by Francisco Simoes, March 1999.

  REAL(KIND=8), INTENT(IN) :: d,nu,s,g,p,csf
  REAL(KIND=8) :: dstar,logdstar,r1,r2,r3,tanhdstar,wstar
  REAL(KIND=8), PARAMETER :: one=1.0

  dstar = (s-one)*g*d**3/(nu*nu)
  logdstar = LOG10(dstar)
  tanhdstar = DTANH(logdstar-4.6D0)
  r1 = -3.76715D0+1.92944*logdstar-9.815D-2*logdstar*logdstar- &
       5.75D-3*logdstar**3+5.6D-4*logdstar**4
  r2 = LOG10(one-(one-csf)/8.5D-1)
  r2 = r2-tanhdstar*(one-csf)**2.3
  r2 = r2+0.3D0*(0.5D0-csf)*(one-csf)*(one-csf)*(logdstar-4.6D0)
  r3 = (6.5D-1-csf*tanhdstar/2.83D0)**(one+(3.5D0-p)/2.5D0)
  wstar = r3*10.0D0**(r1+r2)
  falveld = (wstar*(s-one)*g*nu)**(one/3.0D0)

END FUNCTION falveld
