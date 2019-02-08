MODULE shr_kind_mod

  !----------------------------------------------------------------------------
  ! precision/kind constants add data public
  !----------------------------------------------------------------------------
  use ISO_C_BINDING, only: C_DOUBLE, C_INT, C_FLOAT, C_BOOL
  public
  integer,parameter :: SHR_KIND_R8 = C_DOUBLE ! 8 byte real
  integer,parameter :: SHR_KIND_R4 = C_FLOAT ! 4 byte real
  integer,parameter :: SHR_KIND_RN = kind(1.0)              ! native real
  integer,parameter :: SHR_KIND_I8 = selected_int_kind (13) ! 8 byte integer
  integer,parameter :: SHR_KIND_I4 = C_INT ! 4 byte integer
  integer,parameter :: SHR_KIND_IN = kind(1)                ! native integer
  integer,parameter :: SHR_KIND_CS = 80                     ! short char
  integer,parameter :: SHR_KIND_CM = 160                    ! mid-sized char
  integer,parameter :: SHR_KIND_CL = 256                    ! long char
  integer,parameter :: SHR_KIND_CX = 512                    ! extra-long char
  integer,parameter :: SHR_KIND_CXX= 4096                   ! extra-extra-long char
  integer,parameter :: SHR_KIND_BOOL = C_BOOL               ! boolean
  
END MODULE shr_kind_mod
