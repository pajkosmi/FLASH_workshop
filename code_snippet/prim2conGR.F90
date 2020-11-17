!!****if* source/physics/Hydro/HydroMain/Spark/prim2con
!!
!! NAME
!!
!!  hy_rk_prim2con
!!
!! SYNOPSIS
!!
!!  hy_rk_prim2con( real(IN)  :: V(HY_VARINUM3),
!!                   real(OUT) :: CU(HY_VARINUM))
!!
!! ARGUMENTS
!!
!! V  - primitive variables + GAMC,GAME
!! CU - conservative variables
!! vec_size - size of primitive vector
!! update_Soln - Is this being called from hy_rk_updateSoln?
!! T - stress energy tensor
!!
!! DESCRIPTION
!!
!!  This routine calculates conversions from primitive variables to conservative variables for GR Hydro.  See Toro or Rezzola for a review
!!  GRMHD routine from Gammie+ 2003 & Mosta 2013
!!  Also optionally calculates stress energy tensor
!!
!!***

subroutine prim2con(V,CU,vec_size,update_Soln,T)

  use Hydro_data, ONLY : hy_c2g_dens, hy_c2g_pres, hy_c2g_eint,&
                         hy_c2g_c, hy_c2g_B
  implicit none

#include "Flash.h"
#include "Spark.h"

  !! Arguments type declaration -----------
  integer, intent(in)                       :: vec_size
  !real ,dimension(HY_NUM_VARS), intent(IN)  :: V !OG
  real ,dimension(vec_size), intent(IN)     :: V
  real ,dimension(HY_NUM_FLUX), intent(OUT) :: CU
  logical,             intent(IN), optional :: update_Soln
  real, intent(OUT), optional               :: T(0:3,0:3)
  real                                      :: g(0:3,0:3)
  real                                      :: W, h, hm1, detg, alpha
  real, dimension(3,3)                      :: gamco,gam
  real, dimension(3)                        :: beta, Smag
  real, dimension(0:3)                      :: bc
  !! --------------------------------------
  real  :: u2, bc2, Bv
  real  :: uvel(0:4)
  real  :: f_d, f_P, f_e, f_v, f_B, f_psi !conversion factors
  integer :: vDENS, vPRES, vVELX, vVELY, vVELZ, vMAGX, vMAGY, vMAGZ, vPSIB!indices
  integer :: vG00, vG11, vG22, vG33
  integer :: mu,nu

  !Select indices whether pencils or solnData/hy_starState is being accessed

  !if hy_rk_updateSoln() calls this we need to rename variable indices
  !b/c V is not a pencil
  !Also we need to convert to geometric units (assume the subroutine is #inlined)
  if (present(update_Soln)) then
    !conversion factors
    f_D = hy_c2g_dens
    f_P = hy_c2g_pres
    f_e = hy_c2g_eint
    f_v = hy_c2g_c
    f_B = hy_c2G_B
    f_psi = hy_c2G_B*hy_c2g_c !unsure about factor of sqrt(4PI) 
    !indices
    vDENS = DENS_VAR
    vPRES = PRES_VAR
    vVELX = VELX_VAR
    vVELY = VELY_VAR
    vVELZ = VELZ_VAR
    vMAGX = MAGX_VAR
    vMAGY = MAGY_VAR
    vMAGZ = MAGZ_VAR
    vPSIB = PSIB_VAR
    vG00  = G00_VAR
    vG11  = G11_VAR
    vG22  = G22_VAR
    vG33  = G33_VAR
  else!riemann solver calls this with pencil indices
    !conversion factors (already in geometric)
    f_D = 1.
    f_P = 1.
    f_e = 1.
    f_v = 1.       
    f_B = 1.
    f_psi = 1. 
    !indices
    vDENS = HY_DENS
    vPRES = HY_PRES
    vVELX = HY_VELX
    vVELY = HY_VELY
    vVELZ = HY_VELZ
    vMAGX = HY_MAGX
    vMAGY = HY_MAGY
    vMAGZ = HY_MAGZ
    vPSIB = HY_PSIB
    vG00  = HY_G00
    vG11  = HY_G11
    vG22  = HY_G22
    vG33  = HY_G33
  endif

  CU = 0.0

  !!Metric information
  !Spatial metric (covariant) gamma_ij(flat).  Later this will be input
  gamco = 0; gamco(1,1)=V(vG11); gamco(2,2)=V(vG22); gamco(3,3)=V(vG33)

  !print*,"G11 pencil",V(vG11),V(vDENS)
 
  !" " (contravariant) gamma_ik * gamma^kj = dirac_i^j   
  gam = 0  ; gam(1,1)=1./gamco(1,1) ; gam(2,2)=1./gamco(2,2); gam(3,3)=1./gamco(3,3)
  
  !Square root of determinate of covariant spatial metric.
  !detg = sqrt(gamco(1,1)*gamco(2,2)*gamco(3,3))
  detg = sqrt(gamco(1,1))

  !lapse
  alpha = sqrt(-V(vG00))
  !shift (contravariant) Beta^i
  beta = (/0.,0.,0./)

  !!Velocity information, assume FLASH outputs contravariant components
  u2 = gamco(1,1)*V(vVELX)*V(vVELX) + gamco(2,2)*V(vVELY)*V(vVELY)+&
       gamco(3,3)*V(vVELZ)*V(vVELZ)
  
  !convert to appropriate units if needed 
  u2 = u2*f_v*f_v
  !Lorentzian
  W = (1-u2)**(-0.5)

  !!Magnetic field terms
  bc2  = 0.
  bc   = 0.
  Smag = 0.
#ifdef SPARK_GLM
  !B^i*v_i
  Bv = gamco(1,1)*V(vMAGX)*V(vVELX) + gamco(2,2)*V(vMAGY)*V(vVELY) +&
       gamco(3,3)*V(vMAGZ)*V(vVELZ) 
  
  Bv = Bv*f_B*f_v
  
  !bc^mu = magnetic field in (c)omoving frame (contravariant)
  bc(0)   = W/alpha*Bv
  bc(1:3) = V(vMAGX:vMAGZ)*f_B/W + W*Bv*(V(vVELX:vVELZ)*f_v-beta/alpha)
  !B^i*B_i/W^2 + (B^i*v_i)^2
  bc2 = (gamco(1,1)*V(vMAGX)*V(vMAGX) + gamco(2,2)*V(vMAGY)*V(vMAGY) +&
         gamco(3,3)*V(vMAGZ)*V(vMAGZ))*f_B*f_B/(W*W) + Bv*Bv 

  !Magnetic contribution to relativistic momentum conserved variables
  Smag(1) = -1.*alpha*bc(0)*gamco(1,1)*bc(1)
  Smag(2) = -1.*alpha*bc(0)*gamco(2,2)*bc(2) 
  Smag(3) = -1.*alpha*bc(0)*gamco(3,3)*bc(3)

  !B^i
  CU(HY_FMGX:HY_FMGZ) = detg*V(vMAGX:vMAGZ)*f_B
  !Does this change Mike in GR?
  CU(HY_FPSI) = V(vPSIB)*f_psi
#endif

  !Specific enthalpy = c**2 + specific int. energy + (pressure + mag_pressure)/mass density
  !See Spark/Config

  !      specific int. ener. = (spec. int. ener. * dens)/dens
  if (present(update_Soln)) then
    hm1 = V(EINT_VAR)*f_e + (V(vPRES)*f_P+bc2)/(V(vDENS)*f_D)
  else!Pencils only have RHOE variable (energy density) (already Geo units)
    hm1 = V(HY_RHOE)/V(vDENS) + (V(vPRES) + bc2)/V(vDENS)
  endif

  h = 1 + hm1

  ! *GR change according to Mosta+ 2013
  !D 
  CU(HY_MASS) = detg*W*V(vDENS)*f_D
  !S_j
  CU(HY_XMOM) = detg*(V(vDENS)*h*W*W*V(vVELX)*gamco(1,1)*f_D*f_v + Smag(1))
  CU(HY_YMOM) = detg*(V(vDENS)*h*W*W*V(vVELY)*gamco(2,2)*f_D*f_v + Smag(2))
  CU(HY_ZMOM) = detg*(V(vDENS)*h*W*W*V(vVELZ)*gamco(3,3)*f_D*f_v + Smag(3))
  !tau
  CU(HY_ENER) = CU(HY_MASS)*(W - 1) + CU(HY_MASS)*W*hm1 - &
                detg*(V(vPRES)*f_P + 0.5*bc2 + (alpha*bc(0))**2)

  !!~~
  !Stress energy tensor calculation (optional)
  ![Eqn. (6) of Mosta+ 2013]
  if (present(T)) then
    T = 0.
    !Time component of 4 velocity
    uvel(0) = W/alpha
    !Spatial components of 4 velocity
    uvel(1:3) = W*V(vVELX:vVELZ)*f_v

    g = 0
    !contravariant metric for source terms
    g(0,0) = 1./V(G00_VAR)  
    g(1,1) = gam(1,1)
    g(2,2) = gam(2,2) 
    g(3,3) = gam(3,3) 

    do mu=0,3
      do nu=0,3
        
        T(mu,nu) = V(vDENS)*f_D*h*uvel(mu)*uvel(nu) - bc(mu)*bc(nu) + &
                   (V(vPRES)*f_P + 0.5*bc2)*g(mu,nu)
      enddo
    enddo
  endif
end subroutine prim2con
