module mod_abundtypes

TYPE line
        CHARACTER*11 :: name
        CHARACTER*20 :: ion
        DOUBLE PRECISION :: wavelength
        DOUBLE PRECISION :: intensity
        DOUBLE PRECISION :: int_err
        DOUBLE PRECISION :: freq
        DOUBLE PRECISION :: int_dered
        CHARACTER*20 :: transition
        DOUBLE PRECISION :: abundance
        CHARACTER*4 :: zone !high, medium or low ionisation zone 
        INTEGER :: location !for CELs and He, this indicates the position of the line in the main linelist array, so that when its abundance is calculated it can be copied back into the linelist array for easy outputting
        CHARACTER*15 :: latextext ! ion name in latex format
END TYPE

end module mod_abundtypes

module mod_resultarrays

TYPE resultarray
        double precision :: NC_abund_CEL
        double precision :: cii_abund_CEL
        double precision :: ciii_abund_CEL
        double precision :: civ_abund_CEL
        double precision :: C_icf_CEL
        double precision :: C_abund_CEL
        double precision :: Nii_abund_CEL
        double precision :: Niii_abund_CEL
        double precision :: Niv_abund_CEL
        double precision :: Nv_abund_CEL
        double precision :: N_icf_CEL
        double precision :: N_abund_CEL
        double precision :: NO_abund_CEL
        double precision :: Oii_abund_CEL
        double precision :: Oiii_abund_CEL
        double precision :: Oiv_abund_CEL
        double precision :: O_icf_CEL
        double precision :: O_abund_CEL
        double precision :: Neii_abund_CEL
        double precision :: Neiii_abund_CEL
        double precision :: Neiv_abund_CEL
        double precision :: Nev_abund_CEL
        double precision :: Ne_icf_CEL
        double precision :: Ne_abund_CEL
        double precision :: Ariii_abund_CEL
        double precision :: Ariv_abund_CEL
        double precision :: Arv_abund_CEL
        double precision :: Ar_icf_CEL
        double precision :: Ar_abund_CEL
        double precision :: Sii_abund_CEL
        double precision :: Siii_abund_CEL
        double precision :: Cliii_abund_CEL
        double precision :: Cl_icf_CEL
        double precision :: Cl_abund_CEL
        double precision :: S_icf_CEL
        double precision :: S_abund_CEL
        double precision :: Hei_abund_ORL
        double precision :: Heii_abund_ORL
        double precision :: He_abund_ORL
        double precision :: Cii_abund_ORL
        double precision :: Ciii_abund_ORL
        double precision :: C_icf_ORL
        double precision :: C_abund_ORL
        double precision :: Nii_abund_ORL
        double precision :: Niii_abund_ORL
        double precision :: N_icf_ORL
        double precision :: N_abund_ORL
        double precision :: Oii_abund_ORL
        double precision :: O_icf_ORL
        double precision :: O_abund_ORL
        double precision :: Neii_abund_ORL
        double precision :: Ne_icf_ORL
        double precision :: Ne_abund_ORL
        double precision :: OII_density
        double precision :: SII_density
        double precision :: low_density
        double precision :: OII_temp
        double precision :: NII_temp
        double precision :: SII_temp
        double precision :: OI_temp
        double precision :: CI_temp
        double precision :: low_temp
        double precision :: ClIII_density
        double precision :: ArIV_density
        double precision :: CIII_density
        double precision :: OIII_IR_density
        double precision :: SIII_IR_density
        double precision :: ArIII_IR_density
        double precision :: NeIII_IR_density
        double precision :: med_density
        double precision :: OIII_temp
        double precision :: OIII_IR_temp
        double precision :: OIII_UV_temp
        double precision :: NeIII_temp
        double precision :: NeIII_IR_temp
        double precision :: ArIII_temp
        double precision :: SIII_temp
        double precision :: med_temp
        double precision :: NeIV_density
        double precision :: high_density
        double precision :: ArV_temp
        double precision :: NeV_temp
        double precision :: high_temp
        double precision :: mean_cHb
        double precision :: O_R23_upper
        double precision :: O_R23_lower
        double precision :: O_N2
        double precision :: O_O3N2
        double precision :: O_Ar3O3
        double precision :: O_S3O3
        double precision :: adf_O
        double precision :: adf_O2plus
        double precision :: adf_N
        double precision :: adf_N2plus
        double precision :: adf_C
        double precision :: adf_C2plus
        double precision :: adf_Ne
        double precision :: adf_Ne2plus
        double precision :: OII_density_ratio
        double precision :: SII_density_ratio
        double precision :: NII_temp_ratio
        double precision :: OII_temp_ratio
        double precision :: SII_temp_ratio
        double precision :: CI_temp_ratio
        double precision :: OI_temp_ratio
        double precision :: ClIII_density_ratio
        double precision :: ArIV_density_ratio
        double precision :: CIII_density_ratio
        double precision :: OIII_IR_density_ratio
        double precision :: ArIII_IR_density_ratio
        double precision :: SIII_IR_density_ratio
        double precision :: NeIII_IR_density_ratio
        double precision :: OIII_temp_ratio
        double precision :: NeIII_temp_ratio
        double precision :: ArIII_temp_ratio
        double precision :: SIII_temp_ratio
        double precision :: OIII_IR_temp_ratio
        double precision :: OIII_UV_temp_ratio
        double precision :: NeIII_IR_temp_ratio
        double precision :: NeIV_density_ratio
        double precision :: ArV_temp_ratio
        double precision :: NeV_temp_ratio
        double precision :: Bal_jump_temp
        double precision :: oii_v1_abund_orl
        double precision :: oii_v2_abund_orl
        double precision :: oii_v5_abund_orl
        double precision :: oii_v10_abund_orl
        double precision :: oii_v11_abund_orl
        double precision :: oii_v12_abund_orl
        double precision :: oii_v19_abund_orl
        double precision :: oii_v20_abund_orl
        double precision :: oii_v25_abund_orl
        double precision :: oii_v28_abund_orl
        double precision :: oii_v33_abund_orl
        double precision :: oii_3d4f_abund_orl
        double precision :: nii_v3_abund_orl
        double precision :: nii_v5_abund_orl
        double precision :: nii_v8_abund_orl
        double precision :: nii_v12_abund_orl 
        double precision :: nii_v20_abund_orl 
        double precision :: nii_v28_abund_orl
        double precision :: nii_3d4f_abund_orl
end type

end module mod_resultarrays

module mod_atomicdata

type atomic_data

        integer :: NTEMPS
        integer :: NLEVS
        integer :: irats
        character*10 :: ion
        real*8,allocatable :: temps(:)
        real*8,allocatable :: roott(:)
        character*20,allocatable :: labels(:)
        real*8,allocatable :: G(:)
        real*8,allocatable :: waveno(:)
        real*8,allocatable :: A_coeffs(:,:)
        real*8,allocatable :: col_str(:,:,:)

end type

end module mod_atomicdata
