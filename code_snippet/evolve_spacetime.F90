!!
!! Evolve the spacetime according to Eqns. (4-6) of GR1D paper [O'Connor & Ott 2011]
!!
!!  |        |         |
!!  |        |---dr----|
!!  |        |         |
!!     i-1  L_i   i   R_i
!!
!!
!!

subroutine evolve_spacetime(V,xC,dr,Mtot,Phitot,dt,Phiconst,calc_const)
  use Hydro_data, ONLY : hy_c2g_pres,&
      hy_c2g_dens, hy_c2g_eint, hy_c2g_c, hy_c2g_time

  implicit none

  real :: xC
  real, dimension(NUNK_VARS) :: V
  real :: Mtot, Phitot, dMin, dMout
  real :: h, W, rho, eps, P, velx, frac, dPhi
  real, intent(IN) :: dr, dt
  real, intent(INOUT) :: Phiconst
  logical, intent(IN) :: calc_const

  !Primitives in geo
  P   = V(PRES_VAR)*hy_c2g_pres
  rho = V(DENS_VAR)*hy_c2g_dens
  eps = V(EINT_VAR)*hy_c2g_eint
  velx = V(VELX_VAR)*hy_c2g_c

  !Specific enthalpy
  h = 1 + (eps + P/rho)

  !Lorentz factor
  W = (1 - V(G11_VAR)*velx*velx)**(-0.5)

  !!Cell centered
  !dMint/dr
  V(DMNT_VAR) = (4*PI*xC*xC*(rho*h*W*W - P))

  !Mass_int   
  !start at 0
  !First add left half of cell
  dMin = 4./3*PI*((xC)**3 - (xC-0.5*dr)**3)*(rho*h*W*W - P)
  Mtot = Mtot + dMin
  !Store new enclosed mass value
  V(MINT_VAR) = Mtot
  
  !Next add right half of cell
  dMout = 4./3*PI*((xC+0.5*dr)**3 - (xC)**3)*(rho*h*W*W - P)
  
  Mtot = Mtot + dMout

  !g11 covariant metric term
  V(G11_VAR) = (1. - 2.*V(MINT_VAR)/xC)**(-1)

  !Radial derivative of metric potential
  !dphi/dr Eqn. (A.4) GR1D Paper
  V(DMPT_VAR) = V(G11_VAR)*(V(MINT_VAR)/xC/xC + &
         4*PI*xC*(P + rho*h*W*W*V(G11_VAR)*velx*velx))

  !phi (add previous constant)
  !start at phi0 
  !Add left contribution

  !fractional change in metric potential
  dPhi = V(DMPT_VAR)*dr 

  ! Time derivative of metric potential
  !d/dt(MPOT_VAR) = DPTT_VAR
  V(DPTT_VAR) = (0.5*(2*Phitot + dPhi) - V(MPOT_VAR))/(dt)

  !Metric potential
  !Save phi at cell center by averaging left & right face values
  V(MPOT_VAR) = 0.5*(2*Phitot + dPhi)

  !Now potential at right face
  Phitot = Phitot + dPhi

  !Radial derivative of g11 metric term
  !dg11/dr = 2*X*dXdr = 2*X^4*(dmdr/r - m/r/r)
  V(DG11_VAR) = 2*V(G11_VAR)*V(G11_VAR)*(V(DMNT_VAR)/xC - &
                                        V(MINT_VAR)/xC/xC)

  !Edge cases to take care of integration constant at boundary
  if (calc_const) then
    Phiconst = 0.5*log(1 - 2*Mtot/(xC + 0.5*dr)) - Phitot
  endif
 

  !Right now these are calculated at the bottom of hy_rk_updateSoln()
  !!alpha = exp(phi)
  !V(LAPS_VAR) = exp(V(MPOT_VAR))

  !!g00 = -alph**2
  !V(G00_VAR) = - V(LAPS_VAR)*V(LAPS_VAR)

end subroutine evolve_spacetime
