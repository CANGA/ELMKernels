

void ReadSNICAR(const std::string &file)


subroutine SnowOptics_init( )
     
     use fileutils  , only : getfil
     use elm_varctl , only : fsnowoptics
     use spmdMod    , only : masterproc
     use ncdio_pio  , only : file_desc_t, ncd_io, ncd_pio_openfile, ncd_pio_closefile
     use ncdio_pio  , only : ncd_pio_openfile, ncd_inqfdims, ncd_pio_closefile, ncd_inqdid, ncd_inqdlen

     type(file_desc_t)  :: ncid                        ! netCDF file id
     character(len=256) :: locfn                       ! local filename
     character(len= 32) :: subname = 'SnowOptics_init' ! subroutine name
     integer            :: ier                         ! error status

    !mgf++
    logical :: readvar      ! determine if variable was read from NetCDF file
    !mgf--

     !
     ! Open optics file:
     if(masterproc) write(iulog,*) 'Attempting to read snow optical properties .....'
     call getfil (fsnowoptics, locfn, 0)
     call ncd_pio_openfile(ncid, locfn, 0)
     if(masterproc) write(iulog,*) subname,trim(fsnowoptics)

     ! direct-beam snow Mie parameters:
     call ncd_io('ss_alb_ice_drc', ss_alb_snw_drc,            'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ice_drc',asm_prm_snw_drc,          'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ice_drc', ext_cff_mss_snw_drc, 'read', ncid, posNOTonfile=.true.)

     ! diffuse snow Mie parameters
     call ncd_io( 'ss_alb_ice_dfs', ss_alb_snw_dfs,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ice_dfs', asm_prm_snw_dfs,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ice_dfs', ext_cff_mss_snw_dfs, 'read', ncid, posNOTonfile=.true.)

#ifdef MODAL_AER
    !mgf++
    ! size-dependent BC parameters and BC enhancement factors
    if (masterproc) write(iulog,*) 'Attempting to read optical properties for within-ice BC (modal aerosol treatment) ...'

    ! BC species 1 Mie parameters
    call ncd_io( 'ss_alb_bc_mam', ss_alb_bc1,           'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()
    call ncd_io( 'asm_prm_bc_mam', asm_prm_bc1,         'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()
    call ncd_io( 'ext_cff_mss_bc_mam', ext_cff_mss_bc1, 'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()

    ! BC species 2 Mie parameters (identical, before enhancement factors applied)
    call ncd_io( 'ss_alb_bc_mam', ss_alb_bc2,           'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()
    call ncd_io( 'asm_prm_bc_mam', asm_prm_bc2,         'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()
    call ncd_io( 'ext_cff_mss_bc_mam', ext_cff_mss_bc2, 'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()

    ! size-dependent BC absorption enhancement factors for within-ice BC
    call ncd_io( 'bcint_enh_mam', bcenh, 'read', ncid, readvar=readvar, posNOTonfile=.true.)
    if (.not. readvar) call endrun()

#else
    ! bulk aerosol treatment
     ! BC species 1 Mie parameters
     call ncd_io( 'ss_alb_bcphil', ss_alb_bc1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_bcphil', asm_prm_bc1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_bcphil', ext_cff_mss_bc1, 'read', ncid, posNOTonfile=.true.)

     ! BC species 2 Mie parameters
     call ncd_io( 'ss_alb_bcphob', ss_alb_bc2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_bcphob', asm_prm_bc2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_bcphob', ext_cff_mss_bc2, 'read', ncid, posNOTonfile=.true.)

    !mgf--
#endif

     ! OC species 1 Mie parameters
     call ncd_io( 'ss_alb_ocphil', ss_alb_oc1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ocphil', asm_prm_oc1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ocphil', ext_cff_mss_oc1, 'read', ncid, posNOTonfile=.true.)

     ! OC species 2 Mie parameters
     call ncd_io( 'ss_alb_ocphob', ss_alb_oc2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_ocphob', asm_prm_oc2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_ocphob', ext_cff_mss_oc2, 'read', ncid, posNOTonfile=.true.)

     ! dust species 1 Mie parameters
     call ncd_io( 'ss_alb_dust01', ss_alb_dst1,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust01', asm_prm_dst1,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust01', ext_cff_mss_dst1, 'read', ncid, posNOTonfile=.true.)

     ! dust species 2 Mie parameters
     call ncd_io( 'ss_alb_dust02', ss_alb_dst2,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust02', asm_prm_dst2,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust02', ext_cff_mss_dst2, 'read', ncid, posNOTonfile=.true.)

     ! dust species 3 Mie parameters
     call ncd_io( 'ss_alb_dust03', ss_alb_dst3,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust03', asm_prm_dst3,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust03', ext_cff_mss_dst3, 'read', ncid, posNOTonfile=.true.)

     ! dust species 4 Mie parameters
     call ncd_io( 'ss_alb_dust04', ss_alb_dst4,           'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'asm_prm_dust04', asm_prm_dst4,         'read', ncid, posNOTonfile=.true.)
     call ncd_io( 'ext_cff_mss_dust04', ext_cff_mss_dst4, 'read', ncid, posNOTonfile=.true.)


     call ncd_pio_closefile(ncid)
     if (masterproc) then

        write(iulog,*) 'Successfully read snow optical properties'
        ! print some diagnostics:
        write (iulog,*) 'SNICAR: Mie single scatter albedos for direct-beam ice, rds=100um: ', &
             ss_alb_snw_drc(71,1), ss_alb_snw_drc(71,2), ss_alb_snw_drc(71,3),     &
             ss_alb_snw_drc(71,4), ss_alb_snw_drc(71,5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for diffuse ice, rds=100um: ',     &
             ss_alb_snw_dfs(71,1), ss_alb_snw_dfs(71,2), ss_alb_snw_dfs(71,3),     &
             ss_alb_snw_dfs(71,4), ss_alb_snw_dfs(71,5)
        if (DO_SNO_OC) then
           write (iulog,*) 'SNICAR: Including OC aerosols from snow radiative transfer calculations'
        else
           write (iulog,*) 'SNICAR: Excluding OC aerosols from snow radiative transfer calculations'
        endif

#ifdef MODAL_AER
       !mgf++
       ! unique dimensionality for modal aerosol optical properties
       write (iulog,*) 'SNICAR: Subset of Mie single scatter albedos for BC: ', &
            ss_alb_bc1(1,1), ss_alb_bc1(1,2), ss_alb_bc1(2,1), ss_alb_bc1(5,1), ss_alb_bc1(1,10), ss_alb_bc2(1,10)
       write (iulog,*) 'SNICAR: Subset of Mie mass extinction coefficients for BC: ', &
            ext_cff_mss_bc2(1,1), ext_cff_mss_bc2(1,2), ext_cff_mss_bc2(2,1), ext_cff_mss_bc2(5,1), ext_cff_mss_bc2(1,10),&
            ext_cff_mss_bc1(1,10)
       write (iulog,*) 'SNICAR: Subset of Mie asymmetry parameters for BC: ', &
            asm_prm_bc1(1,1), asm_prm_bc1(1,2), asm_prm_bc1(2,1), asm_prm_bc1(5,1), asm_prm_bc1(1,10), asm_prm_bc2(1,10)
       write (iulog,*) 'SNICAR: Subset of BC absorption enhancement factors: ', &
            bcenh(1,1,1), bcenh(1,2,1), bcenh(1,1,2), bcenh(2,1,1), bcenh(5,10,1), bcenh(5,1,8), bcenh(5,10,8)
       ! test comparison: ncks -H -C -F -d wvl,5 -d ncl_rds,1 -d ice_rds,8 -v ss_alb_bc_mam,asm_prm_bc_mam,ext_cff_mss_bc_mam,bcint_enh_mam snicar_optics_5bnd_mam_c160322.nc
       !mgf--
#else
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic BC: ', &
             ss_alb_bc1(1), ss_alb_bc1(2), ss_alb_bc1(3), ss_alb_bc1(4), ss_alb_bc1(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic BC: ', &
             ss_alb_bc2(1), ss_alb_bc2(2), ss_alb_bc2(3), ss_alb_bc2(4), ss_alb_bc2(5)
#endif

        if (DO_SNO_OC) then
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophillic OC: ', &
                ss_alb_oc1(1), ss_alb_oc1(2), ss_alb_oc1(3), ss_alb_oc1(4), ss_alb_oc1(5)
           write (iulog,*) 'SNICAR: Mie single scatter albedos for hydrophobic OC: ', &
                ss_alb_oc2(1), ss_alb_oc2(2), ss_alb_oc2(3), ss_alb_oc2(4), ss_alb_oc2(5)
        endif
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 1: ', &
             ss_alb_dst1(1), ss_alb_dst1(2), ss_alb_dst1(3), ss_alb_dst1(4), ss_alb_dst1(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 2: ', &
             ss_alb_dst2(1), ss_alb_dst2(2), ss_alb_dst2(3), ss_alb_dst2(4), ss_alb_dst2(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 3: ', &
             ss_alb_dst3(1), ss_alb_dst3(2), ss_alb_dst3(3), ss_alb_dst3(4), ss_alb_dst3(5)
        write (iulog,*) 'SNICAR: Mie single scatter albedos for dust species 4: ', &
             ss_alb_dst4(1), ss_alb_dst4(2), ss_alb_dst4(3), ss_alb_dst4(4), ss_alb_dst4(5)
        write(iulog,*)
     end if

   end subroutine SnowOptics_init
