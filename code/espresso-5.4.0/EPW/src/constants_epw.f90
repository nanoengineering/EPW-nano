  !                                                                            
  ! Copyright (C) 2010-2016 Samuel Ponce', Roxana Margine, Carla Verdi, Feliciano Giustino 
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt . 
  !
  !-----------------------------------------------------------------------
  MODULE constants_epw
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! Mathematical constants
  ! 
  REAL(DP), PARAMETER :: pi     = 3.141592653589793238462643383279502884197169399375105820974944d0
  REAL(DP), PARAMETER :: twopi  = 2.d0 * pi
  REAL(DP), PARAMETER :: fpi    = 4.d0 * pi
  REAL(DP), PARAMETER :: pibytwo=3.141592653589793238462643383279502884197169399375105820974944d0 / 2.d0
  REAL(DP), PARAMETER :: one    = 1.d0
  REAL(DP), PARAMETER :: two    = 2.d0
  REAL(DP), PARAMETER :: zero   = 0.d0
  REAL(DP), PARAMETER :: e2     = 2.0_DP      ! the square of the electron charge
  COMPLEX(DP), PARAMETER :: ci   = (0.d0, 1.d0)
  COMPLEX(DP), PARAMETER :: cone = (1.d0, 0.d0)
  COMPLEX(DP), PARAMETER :: czero = (0.d0, 0.d0)
  !
  ! Unit conversion factors
  !
  REAL(DP), PARAMETER :: bohr     = 0.52917721092d0
  REAL(DP), PARAMETER :: ryd2mev  = 13605.6981d0
  REAL(DP), PARAMETER :: ryd2ev   = 13.6056981d0
  REAL(DP), PARAMETER :: rydcm1   = 13.6056981d0 * 8065.541d0
  REAL(DP), PARAMETER :: bohr2ang = 0.52917721092d0
  REAL(DP), PARAMETER :: ev2cmm1  = 8065.541d0
  REAL(DP), PARAMETER :: kelvin2eV= 8.6173427909d-05
  REAL(DP), PARAMETER :: ryd2ghz  = 3.289828d6
  !
  ! THL define
  REAL(KIND=DP), PARAMETER :: kB      = 6.33362d-6      ! Boltzmann constant [Ry/K]
  REAL(KIND=DP), PARAMETER :: au2Ohm  = 8.2164712d+3    ! AU to SI [Ohm]
  REAL(KIND=DP), PARAMETER :: au2Amp  = 2.3418037d-3    ! AU to SI [Ampere]
  REAL(KIND=DP), PARAMETER :: au2V    = 19.241363d0     ! AU to SI [Volt]
  REAL(KIND=DP), PARAMETER :: au2j    = 2.1798741d-18   ! AU to SI [Jole]
  REAL(KIND=DP), PARAMETER :: au2m    = 5.29177249d-11  ! AU to SI [meter]
  REAL(KIND=DP), PARAMETER :: au2cm   = 5.29177249d-9   ! AU to SI [centimeter]
  REAL(KIND=DP), PARAMETER :: au2nm   = 5.29177249d-2   ! AU to SI [nanometer]
  REAL(KIND=DP), PARAMETER :: au2s    = 4.8377687d-17   ! AU to SI [second]
  REAL(KIND=DP), PARAMETER :: au2fs   = 4.8377687d-2    ! AU to femtosecond
  REAL(KIND=DP), PARAMETER :: ryd2thz = 3289.8280d0     ! Ry to THz
  REAL(KIND=DP), PARAMETER :: au2ps   = 4.8377687d-5    ! AU to picosecond
  !
  END MODULE constants_epw

