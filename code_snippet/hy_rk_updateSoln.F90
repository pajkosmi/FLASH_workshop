!!****if* source/Simulation/SimulationMain/TOVStar (vvv based on below vvv)
!!        source/physics/Hydro/HydroMain/Spark/hy_rk_updateSoln
!!
!!  NAME
!!
!!  hy_rk_updateSoln
!!
!!  SYNOPSIS
!!
!!  call hy_rk_updateSoln ( type(Grid_tile_t) :: blockDesc )
!!
!!  DESCRIPTION
!!  Update solution based on conservative fluxes previously calculated.  Then convert
!!  conservative to primitive variables.
!!
!!  ARGUMENTS
!!  blockDesc-block descriptor
!!
!!***
!!Reorder(4): hy_starState, solnData, hy_fl[xyz]
subroutine hy_rk_updateSoln (blockDesc, dt, dtOld, limits, coeffs)

  use Simulation_data, ONLY : sim_curved_spacetime,sim_NS_radius,sim_metric_pot_constant 
  use Hydro_data, ONLY : hy_threadWithinBlock, hy_starState, &
       hy_smallE, hy_smalldens, hy_geometry,hy_fluxCorrectPerLevel,&
       hy_fluxCorrect, hy_grav, hy_4piGinv, hy_alphaGLM, hy_C_hyp,&
       hy_flx, hy_fly, hy_flz, hy_c2g_time, hy_c2g_length, &
       hy_c2g_dens, hy_c2g_eint, hy_c2g_c, hy_c2g_pres, hy_c2G_B, &
       hy_update_conserved_vars, hy_evolve_spacetime, hy_phi_pot_constant
  use Driver_interface, ONLY : Driver_abortFlash
  use Grid_interface, ONLY : Grid_getCellCoords,Grid_getCellFaceAreas,&
                             Grid_getCellVolumes,Grid_renormAbundance
  use Grid_tile, ONLY : Grid_tile_t
  use hy_rk_interface, ONLY : prim2con
  
  implicit none

#include "Flash.h"
#include "constants.h"
#include "Spark.h"

  type(Grid_tile_t)   :: blockDesc
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
  real, intent(IN) :: dt, dtOld
  real, dimension(3), intent(IN) :: coeffs

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: lo, hi, loGC, hiGC
  integer :: xLoGC,yLoGC,zLoGC,xHiGC,yHiGC,zHiGC 
 
  real, allocatable, dimension(:) :: xCenter, xLeft, xRight, &
                                     yCenter, zCenter

  integer :: i,j,k,n,g

  real, pointer :: solnData(:,:,:,:)
  real, pointer :: V0(:), Vstar(:)

  real, dimension(NFLUXES) :: U0, Ustar !May have to match U(HY_NUM_FLUX) b/c of prim2con()

  real :: dx, dy, dz, del(MDIM)
  real, pointer :: Fm(:), Fp(:), Gm(:), Gp(:), Hm(:), Hp(:)
  real :: flx_weight1, flx_weight2
  real :: eint, ekin, emag
  ! Geometry factors
  real :: facM, facP
  real, dimension(NFLUXES) :: Sgeo, Sgrv, Stot
  real, allocatable, dimension(:,:,:) :: faceAreas
  real, allocatable, dimension(:,:,:) :: cellVolumes
  integer :: isize, jsize, ksize, lev, ind
  real :: dhdt
  logical :: calc_phi_const, update_Soln=.True.
  real :: dens_atm, eint_atm
  real :: Mtot, Phitot

  Mtot = 0
  Phitot = 0

  calc_phi_const = .False.

  blkLimits(:,:)   = blockDesc%limits
  lo(:) = blkLimits(LOW,:)
  hi(:) = blkLimits(HIGH,:)
  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
  loGC(:) = blkLimitsGC(LOW,:)
  hiGC(:) = blkLimitsGC(HIGH,:)  

  !convenience indices
  xLoGC = loGC(IAXIS); xHiGC = hiGC(IAXIS)
  yLoGC = loGC(JAXIS); yHiGC = hiGC(JAXIS)
  zLoGC = loGC(KAXIS); zHiGC = hiGC(KAXIS)

  iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1

  nullify(solnData)  
  call blockDesc%getDataPtr(solnData,CENTER)

  call blockDesc%deltas(del)
  dx = del(IAXIS)*hy_c2g_length; dy = del(JAXIS)*hy_c2g_length; dz = del(KAXIS)*hy_c2g_length
  dhdt = minval(del(1:NDIM)*hy_c2g_length)/(coeffs(3)*dt*hy_c2g_time)

  if (hy_geometry /= CARTESIAN) then
     lev = blockDesc%level
     
     allocate(faceAreas(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
     call Grid_getCellFaceAreas(IAXIS,lev,loGC,hiGC,faceAreas)
     faceAreas = faceAreas*hy_c2g_length*hy_c2g_length
 
     allocate(cellVolumes(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
     call Grid_getCellVolumes(lev,loGC,hiGC,cellVolumes)
     cellVolumes = cellVolumes*hy_c2g_length*hy_c2g_length*hy_c2g_length

     allocate(xCenter(xLoGC:xHiGC))
     allocate(xLeft(xLoGC:xHiGC))
     allocate(xRight(xLoGC:xHiGC))
     allocate(yCenter(yLoGC:yHiGC))
     allocate(zCenter(zLoGC:zHiGC))

     call Grid_getCellCoords(IAXIS, CENTER, lev, loGC, hiGC, xCenter)
     call Grid_getCellCoords(IAXIS, LEFT_EDGE, lev, loGC, hiGC, xLeft)
     call Grid_getCellCoords(IAXIS, RIGHT_EDGE, lev, loGC, hiGC, xRight) 
     call Grid_getCellCoords(JAXIS, CENTER, lev, loGC, hiGC, yCenter)
     call Grid_getCellCoords(KAXIS, CENTER, lev, loGC, hiGC, zCenter)

     xCenter = xCenter*hy_c2g_length
     xRight  = xRight*hy_c2g_length
     xLeft   = xLeft*hy_c2g_length
     !Currently unused
     !yCenter = yCenter*hy_c2g_length
     !zCenter = zCenter*hy_c2g_length

  endif
 
  !$omp parallel if (hy_threadWithinBlock) &
  !$omp default(none) &
  !$omp private(n,i,j,k,Fm,Fp,Gm,Gp,Hm,Hp,V0,Vstar,U0,Ustar,Sgeo,facM,facP,&
  !$omp         Sgrv,Stot,emag,ekin)&
  !$omp shared(dt,solnData,limits,hy_starState,hy_flx,hy_fly,hy_flz,hy_grav,&
  !$omp        dtOld,xCenter,xLeft,xRight,yCenter,zCenter,&
  !$omp        blkLimits,blkLimitsGC,coeffs,dx,dy,dz,hy_alphaGLM, hy_C_hyp,&
  !$omp        dhdt, hy_smalldens, hy_smallE)

  !  Begin loop over zones
  !$omp do schedule(guided) collapse(3)
  do k = limits(LOW,KAXIS), limits(HIGH,KAXIS)
     do j = limits(LOW,JAXIS), limits(HIGH,JAXIS)
        do i = limits(LOW,IAXIS), limits(HIGH,IAXIS)
           ! Point to old/intermediate states
           V0    => solnData(i,j,k,:)
           Vstar => hy_starState(i,j,k,:)

           ! Point to the correct fluxes
           Fm => hy_flx(i  ,j  ,k  ,:)
           Fp => hy_flx(i+1,j  ,k  ,:)
#if NDIM > 1
           Gm => hy_fly(i  ,j  ,k  ,:)
           Gp => hy_fly(i  ,j+1,k  ,:)
#if NDIM == 3
           Hm => hy_flz(i  ,j  ,k  ,:)
           Hp => hy_flz(i  ,j  ,k+1,:)
#endif
#endif           

           ! Construct vectors of conserved variables
           call prim2con(V0,U0,NPROP_VARS,update_soln)
           U0(HY_NUM_FLUX+1:NFLUXES) = V0(SPECIES_BEGIN:MASS_SCALARS_END)*V0(DENS_VAR)!*W?
 
           call prim2con(Vstar,Ustar,NPROP_VARS,update_soln) 
           Ustar(HY_NUM_FLUX+1:NFLUXES) = Vstar(SPECIES_BEGIN:MASS_SCALARS_END)*Vstar(DENS_VAR)

           if (.NOT.sim_curved_spacetime) then
             ! Flat spacetime
             ! Get geometric factors and source
             call geoFacs(i,j,k,facM,facP,Sgeo,Ustar,Vstar)
             
             ! Get gravitational source terms
             call gravSources(Ustar,hy_grav(:,i,j,k),Sgrv)

             ! Sum total source terms
             Stot = Sgeo + Sgrv

           else
             !Curved spacetime
             !Geometric factors to weight fluxes
             call geoFacs(i,j,k,facM,facP,Sgeo,Ustar,Vstar)
             Sgeo = 0
             !!~~
             !Source terms accounting for curvature coupling to matter (Mosta+ 2013)
             call genrelSources(Vstar,Ustar,Stot,update_soln,xCenter(i))
             Stot = Stot + Sgeo
           endif


           !Dont evolve beyond the atmosphere atmosphere
           if (xCenter(i) <= sim_NS_radius*hy_c2g_length) then
             !Help for debugging.  Remove eventually
             if (hy_update_conserved_vars) then
               ! Now update conserved vector with flux gradients
               Ustar = coeffs(1)*U0 + coeffs(2)*Ustar +coeffs(3)*( &
                    -dt*hy_c2g_time/dx*(facP*Fp-facM*Fm) &
#if NDIM > 1
                    -dt*hy_c2g_time/dy*(Gp-Gm) &
#if NDIM ==3
                    -dt*hy_c2g_time/dz*(Hp-Hm) &
#endif
#endif
                    +dt*hy_c2g_time*Stot)
             endif
           endif
           ! Update primitive variables
           emag = 0.0
#ifdef SPARK_GLM
           !Evolve PsiB
           !                                      hy_C_hyp = 1. (speed of light)
           Vstar(PSIB_VAR) = Ustar(HY_FPSI)*exp(-hy_alphaGLM*1./dhdt)
           !Convert to proper units
           Vstar(PSIB_VAR) = Vstar(PSIB_VAR)/(hy_c2g_c*hy_c2G_B)
#endif
           !Conservative variables to primitive transformation
           call con2primGR(Vstar,Ustar,i,j,k)

           if (xCenter(i) > 0) then
             !evolve spacetime in between stages
             !if boundary, recalculate integration constant
             if (i == limits(HIGH,IAXIS)) calc_phi_const = .True.
             !!~~
             call evolve_spacetime(Vstar,xCenter(i),dx,Mtot,Phitot,dt*hy_c2g_time,&
                                   hy_phi_pot_constant,calc_phi_const)
           endif !(xC(i) > 0)

           ! Divide partial densities by new mass densities to finalize
           ! update of new mass fractions.
           Vstar(SPECIES_BEGIN:MASS_SCALARS_END) = Ustar(HY_NUM_FLUX+1:NFLUXES)/Vstar(DENS_VAR)
#ifdef GPOT_VAR
           ! Now extrapolate gravity to intermediate time state
           ! the star state GPOT_VAR will be reset so that GPOL_VAR isn't screwed up
           Vstar(GPOT_VAR) = coeffs(1)*V0(GPOT_VAR) + coeffs(2)*Vstar(GPOT_VAR) &
                + coeffs(3)*dt*(V0(GPOT_VAR) - V0(GPOL_VAR))/dtOld
#endif
           ! Release pointers
           nullify(V0)
           nullify(Vstar)
           nullify(Fm)
           nullify(Fp)
           nullify(Gm)
           nullify(Gp)
           nullify(Hm)
           nullify(Hp)

        enddo !i
     enddo !j
  enddo !k
  !$omp end do nowait
  !$omp end parallel

  !!~~
  !Shift metric potential by integration constant
  hy_starState(:,:,:,MPOT_VAR) = hy_starState(:,:,:,MPOT_VAR) + hy_phi_pot_constant
  !Lapse function
  hy_starState(:,:,:,LAPS_VAR) = exp(hy_starState(:,:,:,MPOT_VAR))
  !G00 covariant metric term
  hy_starState(:,:,:,G00_VAR ) = -1*hy_starState(:,:,:,LAPS_VAR)*hy_starState(:,:,:,LAPS_VAR)

  call blockDesc%releaseDataPtr(solnData,CENTER)

#if NSPECIES>0
  !Properly normalize species after the update
  solnData => hy_starState
  call Grid_renormAbundance(blockDesc,blkLimitsGC,solnData)
  nullify(solnData)
#endif 

  if (hy_geometry /= CARTESIAN) then
     deallocate(xCenter)
     deallocate(xLeft)
     deallocate(xRight)
     deallocate(yCenter)
     deallocate(zCenter)
     deallocate(faceAreas)
     deallocate(cellVolumes)
  end if

contains
  !!Account for multiplicative factors at each cell face to account for different 
  !!geometries.
  subroutine  geoFacs(i,j,k,facM,facP,Sgeo,U,V)
    implicit none
    integer, intent(IN) :: i,j,k
    real, intent(OUT) :: facM, facP
    real, dimension(NFLUXES) :: Sgeo
    real, intent(IN) :: U(:)
    real, pointer, intent(IN) :: V(:)

    real    :: presStar, densStar, pmomStar, tmomStar, xmomStar
    real    :: pmagStar, xmagStar, zmagStar
    integer :: VEL_PHI, MOM_PHI, MOM_PHI_FLUX, MAG_PHI,  MAG_PHI_FLUX
    integer :: VEL_ZI, MOM_ZI, MOM_ZI_FLUX, MAG_ZI,  MAG_ZI_FLUX
    integer :: VEL_THT, MOM_THT, MOM_THT_FLUX
    real    :: alpha, dx_sph

    if (hy_geometry == CARTESIAN) then
       facM = 1.0; facP = 1.0; Sgeo = 0.0
       return
    endif

    select case(hy_geometry) ! First, select whether y or z is phi-direction
    case(CYLINDRICAL)
       MOM_PHI = HY_ZMOM
       MOM_PHI_FLUX = HY_ZMOM
       MOM_ZI       = HY_YMOM
       MOM_ZI_FLUX  = HY_YMOM
#ifdef SPARK_GLM
       MAG_PHI      = HY_MAGZ
#endif
       alpha = 1.
    case(POLAR)
       MOM_PHI      = HY_YMOM
       MOM_PHI_FLUX = HY_YMOM
       MOM_ZI       = HY_ZMOM
       MOM_ZI_FLUX  = HY_ZMOM
#ifdef SPARK_GLM
       MAG_PHI      = HY_MAGY
#endif
       alpha = 1.
    case(SPHERICAL)
       MOM_PHI      = HY_ZMOM
       MOM_PHI_FLUX = HY_ZMOM
       MOM_THT      = HY_YMOM
       MOM_THT_FLUX = HY_YMOM
       dx_sph = (xRight(i)**3 - xLeft(i)**3) / (3.*xCenter(i)**2)
       alpha  = 2.
    end select
    
    facM = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
    facP = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)

    Sgeo = 0.
    !! Calculate geometrical source terms.  See S&O 75.
    !! alpha 2alpha(GR)*P/X/r?
    Sgeo(HY_XMOM) = (V(DENS_VAR)*hy_c2g_dens*V(VELZ_VAR)*V(VELZ_VAR)*hy_c2g_c*hy_c2g_c + &
                    alpha*V(PRES_VAR)*hy_c2g_pres) / xCenter(i)!T phi,phi

    Sgeo(MOM_PHI) = V(DENS_VAR)*hy_c2g_dens*V(VELZ_VAR)*V(VELX_VAR)*hy_c2g_c*hy_c2g_c / &
                    xCenter(i)!T phi,r

#ifdef SPARK_GLM
    !Mike need to rectify MAGP factor with an emag value!!
    ! P* is the total Pressure
    ! This presently does not work for POLAR coordinates
    Sgeo(HY_XMOM) = Sgeo(HY_XMOM) - ((V(MAGZ_VAR)*hy_c2G_B)**2 - alpha*V(MAGP_VAR)*1)/ xCenter(i)
    Sgeo(MOM_PHI) = Sgeo(MOM_PHI) - V(MAGZ_VAR)*V(MAGX_VAR)*hy_c2G_B*hy_c2G_B / xCenter(i)

    !Mike this is a hack MAGPHI isn't defined for spherical case
    if (hy_geometry /= SPHERICAL) then
      Sgeo(MAG_PHI) = - (V(VELZ_VAR)*V(MAGX_VAR) - V(MAGZ_VAR)*V(VELX_VAR)) / xCenter(i) !O phi,r
      Sgeo(MAG_PHI) = Sgeo(MAG_PHI)*hy_c2G_B*hy_c2g_c
    endif
#endif
    Sgeo(MOM_PHI) = - Sgeo(MOM_PHI)

    if (hy_geometry == SPHERICAL) then
       Sgeo(HY_XMOM) = Sgeo(HY_XMOM) + U(MOM_THT)**2/(V(DENS_VAR)*hy_c2g_dens) / xCenter(i)
       Sgeo = Sgeo*dx/dx_sph
#ifdef SPARK_GR 
       !Needed for 1D spherical formulation
       facP = (xRight(i)/xCenter(i))**2
       facM = (xLeft(i)/ xCenter(i))**2
#endif

    endif
    if (xCenter(i) < 0.0) then
       facM = 0.
       facP = 0.
       Sgeo = 0.
    end if
  end subroutine geoFacs

  !!Calculate source terms due to Newtonian gravity.
  subroutine gravSources(U,g,S)
    implicit none
    real, intent(IN) :: U(NFLUXES)
    real, intent(IN) :: g(MDIM)
    real, intent(OUT) :: S(NFLUXES)
    real :: hy_c2g_accel
    S = 0.
#ifdef GRAVITY
    hy_c2g_accel = hy_c2g_c/hy_c2g_time
    S(HY_XMOM:HY_ZMOM) = U(HY_MASS)*g(:)*hy_c2g_accel
    S(HY_ENER) = dot_product(U(HY_XMOM:HY_ZMOM),g(:)*hy_c2g_accel)
#endif
  end subroutine gravSources

  !!Source terms from GR coupling of spacetime to stress-energy tensor
  subroutine genrelSources(V,U,S,update_soln,rad)
    implicit none
    real, pointer, intent(IN) :: V(:)
    real, intent(IN) :: rad
    real, dimension(NUNK_VARS) :: Vtest
    real, dimension(NFLUXES) :: U, Utest
    real, intent(OUT) :: S(NFLUXES)
    real :: T(0:3,0:3), g(0:3,0:3), GAM(0:3,0:3,0:3),dgdx(0:3,0:3,0:3),dlnadx(0:3)
    real :: detg, alph, X, rho, h, W2, velr, dS
    logical :: update_soln
    integer :: mu, nu

    !Initialize
    S      = 0. !Source term
    T      = 0. !Stress energy (SE) tensor
    GAM    = 0. !Christoffel symbol
    dgdx   = 0. !Derivative of covariant metric w.r.t. direction
    dlnadx = 0. !Derivative of natural log of lapse w.r.t. direction

    !!~~
    !Assign values to Christoffel symbols and metric derivatives
    call calc_metric_derivatives(V,U,rad,GAM,dgdx,dlnadx)  

    !!~~
    !calc contravariant Tmunu
    call prim2con(V,U,NPROP_VARS,update_soln,T)  

    !define covariant metric (1D only)
    g = 0
    g(0,0) = V(G00_VAR)
    g(1,1) = V(G11_VAR)
    g(2,2) = V(G22_VAR)
    g(3,3) = V(G33_VAR)

    !calc det(g), determinate of spatial metric
    !detg = sqrt(g(1,1)*g(2,2)*g(3,3))
    detg = sqrt(g(1,1))    
    alph = V(LAPS_VAR) 

    !only for 1D problem
    X = sqrt(g(1,1))

    !calc source terms numerically

    !Xmomentum [Eqn. (17) Mosta+ 2013]
    do mu = 0,3
      do nu = 0,3
        dS = T(mu,nu)*(dgdx(mu,nu,1) - DOT_PRODUCT(GAM(:,mu,nu),g(:,1)))
        S(HY_XMOM) = S(HY_XMOM) + dS
      enddo 
    enddo

    S(HY_XMOM) = alph*detg*S(HY_XMOM)
  
    !Energy 
    do mu = 0,3
      dS = alph*T(mu,0)*dlnadx(mu)
      do nu = 0,3
        dS = dS - alph*T(mu,nu)*GAM(0,mu,nu)
      enddo
      S(HY_ENER) = S(HY_ENER) + dS
    enddo

    S(HY_ENER) = alph*detg*S(HY_ENER)

    !!!Analytic source terms (old for reference)
    !!!Momentum
    !!!X
    !!!Mosta algebra
    !!!S_XMOM      alpha*detg*(-T00*alph*alph*dphidr 
    !!!                       + T11*(dg11dr - g11*g11/r*(dmdr - m/r)) 
    !!!                       + T22*r 
    !!!                       + T33*r*sin^2(theta))

    !!S(HY_XMOM) = alph*detg*(-1*T(0,0)*alph*alph*V(DMPT_VAR) + &
    !!                         T(1,1)*0.5*V(DG11_VAR) + &
    !!                           T(2,2)*rad + &          !sqrt(G22) = r 
    !!                           T(3,3)*g(3,3)/rad)      !G33/r = r*sin^2(theta)

    !!!No Y/Z source terms for now
    !!!Y
    !!!S(HY_YMOM) = ...

    !!!Z

    !!rho = V(DENS_VAR)*hy_c2g_dens
    !!!Specific enthalpy
    !!h   = 1 + (V(EINT_VAR) + V(PRES_VAR)/V(DENS_VAR))*hy_c2g_eint
    !!!radial velocity
    !!velr = V(VELX_VAR)*hy_c2g_c
    !!!Lorentz factor squared
    !!W2 = 1/(1-V(G11_VAR)*velr*velr)    

    !!!Energy
    !!!          = -alpha*(detg^3)*Mint/r^2*rho*h*W^2*v

    !!S(HY_ENER) = -alph*(detg*detg*detg)*V(MINT_VAR)/rad/rad* &
    !!              rho*h*W2*velr

  end subroutine genrelSources
 
  subroutine calc_metric_Derivatives(V,U,rad,GAM,dgdx,dlnadx)
    implicit none
    real, pointer, intent(IN) :: V(:)
    real, intent(IN) :: rad
    real, dimension(NFLUXES) :: U
    real :: g(0:3,0:3), GAM(0:3,0:3,0:3),dgdx(0:3,0:3,0:3),dlnadx(0:3)
    real :: detg, alph, X, rho, h, W2, velr, dmdt, dphidt

    X = sqrt(V(G11_VAR))
    alph = V(LAPS_VAR)

    rho = V(DENS_VAR)*hy_c2g_dens
    !Specific enthalpy
    h   = 1 + (V(EINT_VAR) + V(PRES_VAR)/V(DENS_VAR))*hy_c2g_eint
    !radial velocity
    velr = V(VELX_VAR)*hy_c2g_c
    !Lorentz factor squared
    W2 = 1/(1-V(G11_VAR)*velr*velr) 

    !Derivative of internal mass w.r.t. time
    dmdt = -4*PI*rad*rad*alph*rho*h*W2*velr

    !Derivative of metric pot. w.r.t. time
    dphidt = V(DPTT_VAR)

    !Indexing explained
    !Christoffel symbol: GAM(upper,lower left,lower right) (symmetric in lower indices)
    ![Values from O'Connor & Ott 2011 Table A.1]
    GAM(0,0,0) = dphidt                   !G^t_tt 
    GAM(0,0,1) = V(DMPT_VAR)              !G^t_tr
    GAM(0,1,0) = GAM(0,0,1)
    GAM(0,1,1) = alph**(-2)*X**4/rad*dmdt !G^t_rr
    GAM(1,0,0) = alph*alph/X/X*V(DMPT_VAR)!G^r_tt
    GAM(1,0,1) = X*X/rad*dmdt             !G^r_tr
    GAM(1,1,0) = GAM(1,0,1)
    GAM(1,1,1) = X*X/rad*(V(DMNT_VAR) - V(MINT_VAR)/rad) !G^r_rr
    GAM(1,2,2) = -rad/X/X                 !G^r_\theta\theta
    GAM(1,3,3) = -V(G33_VAR)/rad/X/X       !G^r_\phi\phi
    GAM(2,1,2) = 1/rad                    !G^\theta_r\theta
    GAM(2,2,1) = GAM(2,1,2)
    GAM(2,3,3) = 0. !-sin(theta)cos(theta)!G^\theta_\phi\phi
    GAM(3,1,3) = 1/rad                    !G^\phi_r\phi
    GAM(3,3,1) = GAM(3,1,3)
    GAM(3,2,3) = 0. !cos(theta)/sin(theta)!G^\phi_\theta\phi
    GAM(3,3,2) = GAM(3,2,3)

    !dgdx(deriv direction,metric ind. 1, metric ind. 2)
    dgdx(0,1,1) = 2*X**4/rad*dmdt !dg11/dt = 2X*dXdt = 2X^4/rad*dmdt

    dgdx(1,1,1) = V(DG11_VAR) !dg11/dr

    !dlnadx(deriv. direction)
    dlnadx(0) = dphidt      !dln(alpha)dt = d(phi)dt

    dlnadx(1) = V(DMPT_VAR) !dln(alpha)dr = d(phi)dr

  end subroutine calc_metric_derivatives

end subroutine hy_rk_updateSoln

!New GR Hydro con2prim routine w/ zerofinder
#include "prim2conGR.F90"
#include "evolve_spacetime.F90"
#ifdef SPARK_GLM
#include "con2primGRMHD.F90"
#else
#include "con2primGR.F90"
#endif
