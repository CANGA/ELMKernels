


namespace ELM {


void SNICAR_RT() {

  int numrad_snw = 5;

  // Define constants
  pi = SHR_CONST_PI;

  // always use Delta approximation for snow
  DELTA = 1;

  if (!Land.urbpoi) {

    // Zero absorbed radiative fluxes:
    for (int i = 0; i <= nlevsno; ++i) {
      flx_abs_lcl(:,:)   = 0._r8
      flx_abs(c_idx,i,:) = 0._r8
    }

    // set snow/ice mass to be used for RT:
    h2osno_lcl = h2osno;



    // Qualifier for computing snow RT: 
    //  1) sunlight from atmosphere model 
    //  2) minimum amount of snow on ground. 
    //     Otherwise, set snow albedo to zero
    if (coszen > 0.0 && h2osno_lcl > min_snw) {     

      // Set variables specific to ELM
      // If there is snow, but zero snow layers, we must create a layer locally.
      // This layer is presumed to have the fresh snow effective radius.
      if (snl < 1) {
         flg_nosnl = 1;
         snl_lcl = 1;
         h2osno_ice_lcl(0) =  h2osno_lcl
         h2osno_liq_lcl(0) =  0._r8
         snw_rds_lcl(0)    =  nint(snw_rds_min)
      } else {
         flg_nosnl         =  0
         snl_lcl           =  snl(c_idx)
         h2osno_liq_lcl(:) =  h2osno_liq(c_idx,:)
         h2osno_ice_lcl(:) =  h2osno_ice(c_idx,:)
         snw_rds_lcl(:)    =  snw_rds(c_idx,:)
      }

               snl_btm   = 0
               snl_top   = snl_lcl+1

               ! for debugging only
               l_idx     = col_pp%landunit(c_idx)
               g_idx     = col_pp%gridcell(c_idx)
               sfctype   = lun_pp%itype(l_idx)
               lat_coord = grc_pp%latdeg(g_idx)
               lon_coord = grc_pp%londeg(g_idx)


#ifdef MODAL_AER
          !mgf++
          !
          ! Assume fixed BC effective radii of 100nm. This is close to
          ! the effective radius of 95nm (number median radius of
          ! 40nm) assumed for freshly-emitted BC in MAM.  Future
          ! implementations may prognose the BC effective radius in
          ! snow.
          rds_bcint_lcl(:)  =  100._r8
          rds_bcext_lcl(:)  =  100._r8
          !mgf--
#endif

            ! Set local aerosol array
            do j=1,sno_nbr_aer
               mss_cnc_aer_lcl(:,j) = mss_cnc_aer_in(c_idx,:,j)
            enddo


            ! Set spectral underlying surface albedos to their corresponding VIS or NIR albedos
            albsfc_lcl(1)                       = albsfc(c_idx,1)
            albsfc_lcl(nir_bnd_bgn:nir_bnd_end) = albsfc(c_idx,2)


            ! Error check for snow grain size:
            do i=snl_top,snl_btm,1
               if ((snw_rds_lcl(i) < snw_rds_min_tbl) .or. (snw_rds_lcl(i) > snw_rds_max_tbl)) then
                  write (iulog,*) "SNICAR ERROR: snow grain radius of ", snw_rds_lcl(i), " out of bounds."
                  write (iulog,*) "NSTEP= ", nstep
                  write (iulog,*) "flg_snw_ice= ", flg_snw_ice
                  write (iulog,*) "column: ", c_idx, " level: ", i, " snl(c)= ", snl_lcl
                  write (iulog,*) "lat= ", lat_coord, " lon= ", lon_coord
                  write (iulog,*) "h2osno(c)= ", h2osno_lcl
                  call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
               endif
            enddo

            ! Incident flux weighting parameters
            !  - sum of all VIS bands must equal 1
            !  - sum of all NIR bands must equal 1
            !
            ! Spectral bands (5-band case)
            !  Band 1: 0.3-0.7um (VIS)
            !  Band 2: 0.7-1.0um (NIR)
            !  Band 3: 1.0-1.2um (NIR)
            !  Band 4: 1.2-1.5um (NIR)
            !  Band 5: 1.5-5.0um (NIR)
            !
            ! The following weights are appropriate for surface-incident flux in a mid-latitude winter atmosphere
            !
            ! 3-band weights
            if (numrad_snw==3) then
               ! Direct:
               if (flg_slr_in == 1) then
                  flx_wgt(1) = 1._r8
                  flx_wgt(2) = 0.66628670195247_r8
                  flx_wgt(3) = 0.33371329804753_r8
                  ! Diffuse:
               elseif (flg_slr_in == 2) then
                  flx_wgt(1) = 1._r8
                  flx_wgt(2) = 0.77887652162877_r8
                  flx_wgt(3) = 0.22112347837123_r8
               endif

               ! 5-band weights
            elseif(numrad_snw==5) then
               ! Direct:
               if (flg_slr_in == 1) then
                  flx_wgt(1) = 1._r8
                  flx_wgt(2) = 0.49352158521175_r8
                  flx_wgt(3) = 0.18099494230665_r8
                  flx_wgt(4) = 0.12094898498813_r8
                  flx_wgt(5) = 0.20453448749347_r8
                  ! Diffuse:
               elseif (flg_slr_in == 2) then
                  flx_wgt(1) = 1._r8
                  flx_wgt(2) = 0.58581507618433_r8
                  flx_wgt(3) = 0.20156903770812_r8
                  flx_wgt(4) = 0.10917889346386_r8
                  flx_wgt(5) = 0.10343699264369_r8
               endif
            endif

            ! Loop over snow spectral bands
            do bnd_idx = 1,numrad_snw

               mu_not    = coszen(c_idx)  ! must set here, because of error handling
               flg_dover = 1              ! default is to redo
               err_idx   = 0              ! number of times through loop

               do while (flg_dover > 0)

                  ! DEFAULT APPROXIMATIONS:
                  !  VIS:       Delta-Eddington
                  !  NIR (all): Delta-Hemispheric Mean
                  !  WARNING:   DO NOT USE DELTA-EDDINGTON FOR NIR DIFFUSE - this sometimes results in negative albedo
                  !  
                  ! ERROR CONDITIONS:
                  !  Conditions which cause "trip", resulting in redo of RT approximation:
                  !   1. negative absorbed flux
                  !   2. total absorbed flux greater than incident flux
                  !   3. negative albedo
                  !   NOTE: These errors have only been encountered in spectral bands 4 and 5
                  !
                  ! ERROR HANDLING
                  !  1st error (flg_dover=2): switch approximation (Edd->HM or HM->Edd)
                  !  2nd error (flg_dover=3): change zenith angle by 0.02 (this happens about 1 in 10^6 cases)
                  !  3rd error (flg_dover=4): switch approximation with new zenith
                  !  Subsequent errors: repeatedly change zenith and approximations...

                  if (bnd_idx == 1) then
                     if (flg_dover == 2) then
                        APRX_TYP = 3
                     elseif (flg_dover == 3) then
                        APRX_TYP = 1
                        if (coszen(c_idx) > 0.5_r8) then
                           mu_not = mu_not - 0.02_r8
                        else
                           mu_not = mu_not + 0.02_r8
                        endif
                     elseif (flg_dover == 4) then
                        APRX_TYP = 3
                     else
                        APRX_TYP = 1
                     endif

                  else
                     if (flg_dover == 2) then
                        APRX_TYP = 1
                     elseif (flg_dover == 3) then
                        APRX_TYP = 3
                        if (coszen(c_idx) > 0.5_r8) then
                           mu_not = mu_not - 0.02_r8
                        else
                           mu_not = mu_not + 0.02_r8
                        endif
                     elseif (flg_dover == 4) then
                        APRX_TYP = 1
                     else
                        APRX_TYP = 3
                     endif

                  endif

                  ! Set direct or diffuse incident irradiance to 1
                  ! (This has to be within the bnd loop because mu_not is adjusted in rare cases)
                  if (flg_slr_in == 1) then
                     flx_slrd_lcl(bnd_idx) = 1._r8/(mu_not*pi) ! this corresponds to incident irradiance of 1.0
                     flx_slri_lcl(bnd_idx) = 0._r8
                  else
                     flx_slrd_lcl(bnd_idx) = 0._r8
                     flx_slri_lcl(bnd_idx) = 1._r8
                  endif

                  ! Pre-emptive error handling: aerosols can reap havoc on these absorptive bands.
                  ! Since extremely high soot concentrations have a negligible effect on these bands, zero them.
                  if ( (numrad_snw == 5).and.((bnd_idx == 5).or.(bnd_idx == 4)) ) then
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif

                  if ( (numrad_snw == 3).and.(bnd_idx == 3) ) then
                     mss_cnc_aer_lcl(:,:) = 0._r8
                  endif

                  ! Define local Mie parameters based on snow grain size and aerosol species,
                  !  retrieved from a lookup table.
                  if (flg_slr_in == 1) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (direct radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_drc(rds_idx,bnd_idx)
                        asm_prm_snw_lcl(i)     = asm_prm_snw_drc(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_drc(rds_idx,bnd_idx)
                     enddo
                  elseif (flg_slr_in == 2) then
                     do i=snl_top,snl_btm,1
                        rds_idx = snw_rds_lcl(i) - snw_rds_min_tbl + 1
                        ! snow optical properties (diffuse radiation)
                        ss_alb_snw_lcl(i)      = ss_alb_snw_dfs(rds_idx,bnd_idx)
                        asm_prm_snw_lcl(i)     = asm_prm_snw_dfs(rds_idx,bnd_idx)
                        ext_cff_mss_snw_lcl(i) = ext_cff_mss_snw_dfs(rds_idx,bnd_idx)
                     enddo
                  endif

!H. Wang
                  ! aerosol species 1 optical properties
                 ! ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)      
                 ! asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                 ! ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                  ! aerosol species 2 optical properties
                 ! ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)      
                 ! asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                 ! ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
!H. Wang
                  ! aerosol species 3 optical properties
                  ss_alb_aer_lcl(3)        = ss_alb_oc1(bnd_idx)      
                  asm_prm_aer_lcl(3)       = asm_prm_oc1(bnd_idx)
                  ext_cff_mss_aer_lcl(3)   = ext_cff_mss_oc1(bnd_idx)

                  ! aerosol species 4 optical properties
                  ss_alb_aer_lcl(4)        = ss_alb_oc2(bnd_idx)      
                  asm_prm_aer_lcl(4)       = asm_prm_oc2(bnd_idx)
                  ext_cff_mss_aer_lcl(4)   = ext_cff_mss_oc2(bnd_idx)

                  ! aerosol species 5 optical properties
                  ss_alb_aer_lcl(5)        = ss_alb_dst1(bnd_idx)      
                  asm_prm_aer_lcl(5)       = asm_prm_dst1(bnd_idx)
                  ext_cff_mss_aer_lcl(5)   = ext_cff_mss_dst1(bnd_idx)

                  ! aerosol species 6 optical properties
                  ss_alb_aer_lcl(6)        = ss_alb_dst2(bnd_idx)      
                  asm_prm_aer_lcl(6)       = asm_prm_dst2(bnd_idx)
                  ext_cff_mss_aer_lcl(6)   = ext_cff_mss_dst2(bnd_idx)

                  ! aerosol species 7 optical properties
                  ss_alb_aer_lcl(7)        = ss_alb_dst3(bnd_idx)      
                  asm_prm_aer_lcl(7)       = asm_prm_dst3(bnd_idx)
                  ext_cff_mss_aer_lcl(7)   = ext_cff_mss_dst3(bnd_idx)

                  ! aerosol species 8 optical properties
                  ss_alb_aer_lcl(8)        = ss_alb_dst4(bnd_idx)      
                  asm_prm_aer_lcl(8)       = asm_prm_dst4(bnd_idx)
                  ext_cff_mss_aer_lcl(8)   = ext_cff_mss_dst4(bnd_idx)


                  ! 1. snow and aerosol layer column mass (L_snw, L_aer [kg/m^2])
                  ! 2. optical Depths (tau_snw, tau_aer)
                  ! 3. weighted Mie properties (tau, omega, g)

                  ! Weighted Mie parameters of each layer
                  do i=snl_top,snl_btm,1
#ifdef MODAL_AER
                   !mgf++ within-ice and external BC optical properties
                   !
                   ! Lookup table indices for BC optical properties,
                   ! dependent on snow grain size and BC particle
                   ! size.

                   ! valid for 25 < snw_rds < 1625 um:
                   if (snw_rds_lcl(i) < 125) then
                      tmp1 = snw_rds_lcl(i)/50
                      idx_bcint_icerds = nint(tmp1)
                   elseif (snw_rds_lcl(i) < 175) then
                      idx_bcint_icerds = 2
                   else
                      tmp1 = (snw_rds_lcl(i)/250)+2
                      idx_bcint_icerds = nint(tmp1)
                   endif

                   ! valid for 25 < bc_rds < 525 nm
                   idx_bcint_nclrds = nint(rds_bcint_lcl(i)/50)
                   idx_bcext_nclrds = nint(rds_bcext_lcl(i)/50)

                   ! check bounds:
                   if (idx_bcint_icerds < idx_bcint_icerds_min) idx_bcint_icerds = idx_bcint_icerds_min
                   if (idx_bcint_icerds > idx_bcint_icerds_max) idx_bcint_icerds = idx_bcint_icerds_max
                   if (idx_bcint_nclrds < idx_bc_nclrds_min) idx_bcint_nclrds = idx_bc_nclrds_min
                   if (idx_bcint_nclrds > idx_bc_nclrds_max) idx_bcint_nclrds = idx_bc_nclrds_max
                   if (idx_bcext_nclrds < idx_bc_nclrds_min) idx_bcext_nclrds = idx_bc_nclrds_min
                   if (idx_bcext_nclrds > idx_bc_nclrds_max) idx_bcext_nclrds = idx_bc_nclrds_max

                   ! print ice index (debug):
                   !write(iulog,*) "MGF: ice index= ", idx_bcint_icerds

                   ! retrieve absorption enhancement factor for within-ice BC
                   enh_fct = bcenh(bnd_idx,idx_bcint_nclrds,idx_bcint_icerds)

                   ! get BC optical properties (moved from above)
                   ! aerosol species 1 optical properties (within-ice BC)
                   ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx,idx_bcint_nclrds)
                   asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx,idx_bcint_nclrds)
                   ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx,idx_bcint_nclrds)*enh_fct

                   ! aerosol species 2 optical properties (external BC)
                   ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx,idx_bcext_nclrds)
                   asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx,idx_bcext_nclrds)
                   ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx,idx_bcext_nclrds)

#else
                   ! bulk aerosol treatment (BC optical properties independent
                   ! of BC and ice grain size)
                   ! aerosol species 1 optical properties (within-ice BC)
                   ss_alb_aer_lcl(1)        = ss_alb_bc1(bnd_idx)
                   asm_prm_aer_lcl(1)       = asm_prm_bc1(bnd_idx)
                   ext_cff_mss_aer_lcl(1)   = ext_cff_mss_bc1(bnd_idx)

                   ! aerosol species 2 optical properties
                   ss_alb_aer_lcl(2)        = ss_alb_bc2(bnd_idx)
                   asm_prm_aer_lcl(2)       = asm_prm_bc2(bnd_idx)
                   ext_cff_mss_aer_lcl(2)   = ext_cff_mss_bc2(bnd_idx)
#endif
                   !mgf--



                     L_snw(i)   = h2osno_ice_lcl(i)+h2osno_liq_lcl(i)
                     tau_snw(i) = L_snw(i)*ext_cff_mss_snw_lcl(i)

                     do j=1,sno_nbr_aer
                        L_aer(i,j)   = L_snw(i)*mss_cnc_aer_lcl(i,j)
                        tau_aer(i,j) = L_aer(i,j)*ext_cff_mss_aer_lcl(j)
                     enddo

                     tau_sum   = 0._r8
                     omega_sum = 0._r8
                     g_sum     = 0._r8

                     do j=1,sno_nbr_aer
                        tau_sum    = tau_sum + tau_aer(i,j) 
                        omega_sum  = omega_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j))
                        g_sum      = g_sum + (tau_aer(i,j)*ss_alb_aer_lcl(j)*asm_prm_aer_lcl(j))
                     enddo

                     tau(i)    = tau_sum + tau_snw(i)
                     omega(i)  = (1/tau(i))*(omega_sum+(ss_alb_snw_lcl(i)*tau_snw(i)))
                     g(i)      = (1/(tau(i)*omega(i)))*(g_sum+ (asm_prm_snw_lcl(i)*ss_alb_snw_lcl(i)*tau_snw(i)))
                  enddo

                  ! DELTA transformations, if requested
                  if (DELTA == 1) then
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)/(1+g(i))
                        omega_star(i) = ((1-(g(i)**2))*omega(i)) / (1-(omega(i)*(g(i)**2)))
                        tau_star(i)   = (1-(omega(i)*(g(i)**2)))*tau(i)
                     enddo
                  else
                     do i=snl_top,snl_btm,1
                        g_star(i)     = g(i)
                        omega_star(i) = omega(i)
                        tau_star(i)   = tau(i)
                     enddo
                  endif

                  ! Total column optical depth:
                  ! tau_elm(i) = total optical depth above the bottom of layer i
                  tau_elm(snl_top) = 0._r8
                  do i=snl_top+1,snl_btm,1
                     tau_elm(i) = tau_elm(i-1)+tau_star(i-1)
                  enddo

                  ! Direct radiation at bottom of snowpack:
                  F_direct_btm = albsfc_lcl(bnd_idx)*mu_not * &
                       exp(-(tau_elm(snl_btm)+tau_star(snl_btm))/mu_not)*pi*flx_slrd_lcl(bnd_idx)

                  ! Intermediates
                  ! Gamma values are approximation-specific.

                  ! Eddington
                  if (APRX_TYP==1) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = (7-(omega_star(i)*(4+(3*g_star(i)))))/4
                        gamma2(i) = -(1-(omega_star(i)*(4-(3*g_star(i)))))/4
                        gamma3(i) = (2-(3*g_star(i)*mu_not))/4
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 0.5
                     enddo

                     ! Quadrature
                  elseif (APRX_TYP==2) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = (3**0.5)*(2-(omega_star(i)*(1+g_star(i))))/2
                        gamma2(i) = omega_star(i)*(3**0.5)*(1-g_star(i))/2
                        gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 1/(3**0.5)
                     enddo

                     ! Hemispheric Mean
                  elseif (APRX_TYP==3) then
                     do i=snl_top,snl_btm,1
                        gamma1(i) = 2 - (omega_star(i)*(1+g_star(i)))
                        gamma2(i) = omega_star(i)*(1-g_star(i))
                        gamma3(i) = (1-((3**0.5)*g_star(i)*mu_not))/2
                        gamma4(i) = 1-gamma3(i)
                        mu_one    = 0.5
                     enddo
                  endif

                  ! Intermediates for tri-diagonal solution
                  do i=snl_top,snl_btm,1
                     lambda(i) = sqrt(abs((gamma1(i)**2) - (gamma2(i)**2)))
                     GAMMA(i)  = gamma2(i)/(gamma1(i)+lambda(i))

                     e1(i)     = 1+(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                     e2(i)     = 1-(GAMMA(i)*exp(-lambda(i)*tau_star(i)))
                     e3(i)     = GAMMA(i) + exp(-lambda(i)*tau_star(i))
                     e4(i)     = GAMMA(i) - exp(-lambda(i)*tau_star(i))
                  enddo !enddo over snow layers


                  ! Intermediates for tri-diagonal solution
                  do i=snl_top,snl_btm,1
                     if (flg_slr_in == 1) then

                        C_pls_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-(tau_elm(i)+tau_star(i))/mu_not)*   &
                             (((gamma1(i)-(1/mu_not))*gamma3(i))+     &
                             (gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_mns_btm(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-(tau_elm(i)+tau_star(i))/mu_not)*   &
                             (((gamma1(i)+(1/mu_not))*gamma4(i))+     &
                             (gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_pls_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-tau_elm(i)/mu_not)*(((gamma1(i)-(1/mu_not))* &
                             gamma3(i))+(gamma4(i)*gamma2(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                        C_mns_top(i) = (omega_star(i)*pi*flx_slrd_lcl(bnd_idx)* &
                             exp(-tau_elm(i)/mu_not)*(((gamma1(i)+(1/mu_not))* &
                             gamma4(i))+(gamma2(i)*gamma3(i))))/((lambda(i)**2)-(1/(mu_not**2)))

                     else
                        C_pls_btm(i) = 0._r8
                        C_mns_btm(i) = 0._r8
                        C_pls_top(i) = 0._r8
                        C_mns_top(i) = 0._r8
                     endif
                  enddo

                  ! Coefficients for tridiaganol matrix solution
                  do i=2*snl_lcl+1,0,1

                     !Boundary values for i=1 and i=2*snl_lcl, specifics for i=odd and i=even    
                     if (i==(2*snl_lcl+1)) then
                        A(i) = 0
                        B(i) = e1(snl_top)
                        D(i) = -e2(snl_top)
                        E(i) = flx_slri_lcl(bnd_idx)-C_mns_top(snl_top)

                     elseif(i==0) then
                        A(i) = e1(snl_btm)-(albsfc_lcl(bnd_idx)*e3(snl_btm))
                        B(i) = e2(snl_btm)-(albsfc_lcl(bnd_idx)*e4(snl_btm))
                        D(i) = 0
                        E(i) = F_direct_btm-C_pls_btm(snl_btm)+(albsfc_lcl(bnd_idx)*C_mns_btm(snl_btm))

                     elseif(mod(i,2)==-1) then   ! If odd and i>=3 (n=1 for i=3)
                        n=floor(i/2.0)
                        A(i) = (e2(n)*e3(n))-(e4(n)*e1(n))
                        B(i) = (e1(n)*e1(n+1))-(e3(n)*e3(n+1))
                        D(i) = (e3(n)*e4(n+1))-(e1(n)*e2(n+1))
                        E(i) = (e3(n)*(C_pls_top(n+1)-C_pls_btm(n)))+(e1(n)*(C_mns_btm(n)-C_mns_top(n+1)))

                     elseif(mod(i,2)==0) then    ! If even and i<=2*snl_lcl
                        n=(i/2)
                        A(i) = (e2(n+1)*e1(n))-(e3(n)*e4(n+1))
                        B(i) = (e2(n)*e2(n+1))-(e4(n)*e4(n+1))
                        D(i) = (e1(n+1)*e4(n+1))-(e2(n+1)*e3(n+1))
                        E(i) = (e2(n+1)*(C_pls_top(n+1)-C_pls_btm(n)))+(e4(n+1)*(C_mns_top(n+1)-C_mns_btm(n))) 
                     endif
                  enddo

                  AS(0) = A(0)/B(0)
                  DS(0) = E(0)/B(0)

                  do i=-1,(2*snl_lcl+1),-1
                     X(i)  = 1/(B(i)-(D(i)*AS(i+1)))
                     AS(i) = A(i)*X(i)
                     DS(i) = (E(i)-(D(i)*DS(i+1)))*X(i)
                  enddo

                  Y(2*snl_lcl+1) = DS(2*snl_lcl+1)
                  do i=(2*snl_lcl+2),0,1
                     Y(i) = DS(i)-(AS(i)*Y(i-1))
                  enddo

                  ! Downward direct-beam and net flux (F_net) at the base of each layer:
                  do i=snl_top,snl_btm,1
                     F_direct(i) = mu_not*pi*flx_slrd_lcl(bnd_idx)*exp(-(tau_elm(i)+tau_star(i))/mu_not)
                     F_net(i)    = (Y(2*i-1)*(e1(i)-e3(i))) + (Y(2*i)*(e2(i)-e4(i))) + &
                          C_pls_btm(i) - C_mns_btm(i) - F_direct(i)
                  enddo

                  ! Upward flux at snowpack top:
                  F_sfc_pls = (Y(2*snl_lcl+1)*(exp(-lambda(snl_top)*tau_star(snl_top))+ &
                       GAMMA(snl_top))) + (Y(2*snl_lcl+2)*(exp(-lambda(snl_top)* &
                       tau_star(snl_top))-GAMMA(snl_top))) + C_pls_top(snl_top)

                  ! Net flux at bottom = absorbed radiation by underlying surface:
                  F_btm_net = -F_net(snl_btm)


                  ! Bulk column albedo and surface net flux
                  albedo    = F_sfc_pls/((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))
                  F_sfc_net = F_sfc_pls - ((mu_not*pi*flx_slrd_lcl(bnd_idx))+flx_slri_lcl(bnd_idx))

                  trip = 0
                  ! Absorbed flux in each layer
                  do i=snl_top,snl_btm,1
                     if(i==snl_top) then
                        F_abs(i) = F_net(i)-F_sfc_net
                     else
                        F_abs(i) = F_net(i)-F_net(i-1)
                     endif
                     flx_abs_lcl(i,bnd_idx) = F_abs(i)


                     ! ERROR check: negative absorption
                     if (flx_abs_lcl(i,bnd_idx) < -0.00001) then
                        trip = 1
                     endif
                  enddo

                  flx_abs_lcl(1,bnd_idx) = F_btm_net

                  if (flg_nosnl == 1) then
                     ! If there are no snow layers (but still snow), all absorbed energy must be in top soil layer
                     !flx_abs_lcl(:,bnd_idx) = 0._r8
                     !flx_abs_lcl(1,bnd_idx) = F_abs(0) + F_btm_net

                     ! changed on 20070408:
                     ! OK to put absorbed energy in the fictitous snow layer because routine SurfaceRadiation
                     ! handles the case of no snow layers. Then, if a snow layer is addded between now and
                     ! SurfaceRadiation (called in CanopyHydrology), absorbed energy will be properly distributed.
                     flx_abs_lcl(0,bnd_idx) = F_abs(0)
                     flx_abs_lcl(1,bnd_idx) = F_btm_net

                  endif

                  !Underflow check (we've already tripped the error condition above)
                  do i=snl_top,1,1
                     if (flx_abs_lcl(i,bnd_idx) < 0._r8) then
                        flx_abs_lcl(i,bnd_idx) = 0._r8
                     endif
                  enddo

                  F_abs_sum = 0._r8
                  do i=snl_top,snl_btm,1
                     F_abs_sum = F_abs_sum + F_abs(i)
                  enddo


                  !ERROR check: absorption greater than incident flux
                  ! (should make condition more generic than "1._r8")
                  if (F_abs_sum > 1._r8) then
                     trip = 1
                  endif

                  !ERROR check:
                  if ((albedo < 0._r8).and.(trip==0)) then
                     trip = 1
                  endif

                  ! Set conditions for redoing RT calculation 
                  if ((trip == 1).and.(flg_dover == 1)) then
                     flg_dover = 2
                  elseif ((trip == 1).and.(flg_dover == 2)) then
                     flg_dover = 3
                  elseif ((trip == 1).and.(flg_dover == 3)) then
                     flg_dover = 4
                  elseif((trip == 1).and.(flg_dover == 4).and.(err_idx < 20)) then
                     flg_dover = 3
                     err_idx = err_idx + 1
                  elseif((trip == 1).and.(flg_dover == 4).and.(err_idx >= 20)) then
                     flg_dover = 0
                     write(iulog,*) "SNICAR ERROR: FOUND A WORMHOLE. STUCK IN INFINITE LOOP! Called from: ", flg_snw_ice
                     write(iulog,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)
                     write(iulog,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)
                     write(iulog,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                     write(iulog,*) "SNICAR STATS: soot1(0)= ", mss_cnc_aer_lcl(0,1)
                     write(iulog,*) "SNICAR STATS: soot2(0)= ", mss_cnc_aer_lcl(0,2)
                     write(iulog,*) "SNICAR STATS: dust1(0)= ", mss_cnc_aer_lcl(0,3)
                     write(iulog,*) "SNICAR STATS: dust2(0)= ", mss_cnc_aer_lcl(0,4)
                     write(iulog,*) "SNICAR STATS: dust3(0)= ", mss_cnc_aer_lcl(0,5)
                     write(iulog,*) "SNICAR STATS: dust4(0)= ", mss_cnc_aer_lcl(0,6)
                     l_idx     = col_pp%landunit(c_idx)
                     write(iulog,*) "column index: ", c_idx
                     write(iulog,*) "landunit type", lun_pp%itype(l_idx)
                     write(iulog,*) "frac_sno: ", frac_sno(c_idx)
                     call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
                  else
                     flg_dover = 0
                  endif

               enddo !enddo while (flg_dover > 0)

               ! Energy conservation check:
               ! Incident direct+diffuse radiation equals (absorbed+bulk_transmitted+bulk_reflected)
               energy_sum = (mu_not*pi*flx_slrd_lcl(bnd_idx)) + flx_slri_lcl(bnd_idx) - (F_abs_sum + F_btm_net + F_sfc_pls)
               if (abs(energy_sum) > 0.00001_r8) then
                  write (iulog,"(a,e13.6,a,i6,a,i6)") "SNICAR ERROR: Energy conservation error of : ", energy_sum, &
                       " at timestep: ", nstep, " at column: ", c_idx
                  call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
               endif

               albout_lcl(bnd_idx) = albedo

               ! Check that albedo is less than 1
               if (albout_lcl(bnd_idx) > 1.0) then

                  write (iulog,*) "SNICAR ERROR: Albedo > 1.0 at c: ", c_idx, " NSTEP= ",nstep
                  write (iulog,*) "SNICAR STATS: bnd_idx= ",bnd_idx
                  write (iulog,*) "SNICAR STATS: albout_lcl(bnd)= ",albout_lcl(bnd_idx), &
                       " albsfc_lcl(bnd_idx)= ",albsfc_lcl(bnd_idx)
                  write (iulog,*) "SNICAR STATS: landtype= ", sfctype
                  write (iulog,*) "SNICAR STATS: h2osno= ", h2osno_lcl, " snl= ", snl_lcl
                  write (iulog,*) "SNICAR STATS: coszen= ", coszen(c_idx), " flg_slr= ", flg_slr_in

                  write (iulog,*) "SNICAR STATS: soot(-4)= ", mss_cnc_aer_lcl(-4,1)
                  write (iulog,*) "SNICAR STATS: soot(-3)= ", mss_cnc_aer_lcl(-3,1)
                  write (iulog,*) "SNICAR STATS: soot(-2)= ", mss_cnc_aer_lcl(-2,1)
                  write (iulog,*) "SNICAR STATS: soot(-1)= ", mss_cnc_aer_lcl(-1,1)
                  write (iulog,*) "SNICAR STATS: soot(0)= ", mss_cnc_aer_lcl(0,1)

                  write (iulog,*) "SNICAR STATS: L_snw(-4)= ", L_snw(-4)
                  write (iulog,*) "SNICAR STATS: L_snw(-3)= ", L_snw(-3)
                  write (iulog,*) "SNICAR STATS: L_snw(-2)= ", L_snw(-2)
                  write (iulog,*) "SNICAR STATS: L_snw(-1)= ", L_snw(-1)
                  write (iulog,*) "SNICAR STATS: L_snw(0)= ", L_snw(0)

                  write (iulog,*) "SNICAR STATS: snw_rds(-4)= ", snw_rds(c_idx,-4)
                  write (iulog,*) "SNICAR STATS: snw_rds(-3)= ", snw_rds(c_idx,-3)
                  write (iulog,*) "SNICAR STATS: snw_rds(-2)= ", snw_rds(c_idx,-2)
                  write (iulog,*) "SNICAR STATS: snw_rds(-1)= ", snw_rds(c_idx,-1)
                  write (iulog,*) "SNICAR STATS: snw_rds(0)= ", snw_rds(c_idx,0)

                  call endrun(decomp_index=c_idx, elmlevel=namec, msg=errmsg(__FILE__, __LINE__))
               endif

            enddo   ! loop over wvl bands


            ! Weight output NIR albedo appropriately
            albout(c_idx,1) = albout_lcl(1)
            flx_sum         = 0._r8
            do bnd_idx= nir_bnd_bgn,nir_bnd_end
               flx_sum = flx_sum + flx_wgt(bnd_idx)*albout_lcl(bnd_idx)
            end do
            albout(c_idx,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))

            ! Weight output NIR absorbed layer fluxes (flx_abs) appropriately
            flx_abs(c_idx,:,1) = flx_abs_lcl(:,1)
            do i=snl_top,1,1
               flx_sum = 0._r8
               do bnd_idx= nir_bnd_bgn,nir_bnd_end
                  flx_sum = flx_sum + flx_wgt(bnd_idx)*flx_abs_lcl(i,bnd_idx)
               enddo
               flx_abs(c_idx,i,2) = flx_sum / sum(flx_wgt(nir_bnd_bgn:nir_bnd_end))          
            end do

            ! If snow < minimum_snow, but > 0, and there is sun, set albedo to underlying surface albedo
         } elseif ( (coszen(c_idx) > 0._r8) .and. (h2osno_lcl < min_snw) .and. (h2osno_lcl > 0._r8) ) then
            albout(c_idx,1) = albsfc(c_idx,1)
            albout(c_idx,2) = albsfc(c_idx,2)

            ! There is either zero snow, or no sun
         else
            albout(c_idx,1) = 0._r8
            albout(c_idx,2) = 0._r8
         endif    ! if column has snow and coszen > 0

      } enddo    ! loop over all columns

    end associate

  end subroutine SNICAR_RT






}




}