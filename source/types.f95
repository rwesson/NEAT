module mod_abundtypes

TYPE line
        CHARACTER(len=11) :: name
        CHARACTER(len=20) :: ion
        DOUBLE PRECISION :: wavelength
        DOUBLE PRECISION :: wavelength_observed !for use with output from ALFA
        DOUBLE PRECISION :: intensity
        DOUBLE PRECISION :: int_err
        DOUBLE PRECISION :: freq
        DOUBLE PRECISION :: int_dered
        CHARACTER(len=20) :: transition
        DOUBLE PRECISION :: abundance
        CHARACTER(len=4) :: zone !high, medium or low ionisation zone
        INTEGER :: location !for CELs and He, this indicates the position of the line in the main linelist array, so that when its abundance is calculated it can be copied back into the linelist array for easy outputting
        CHARACTER(len=15) :: latextext ! ion name in latex format
        CHARACTER(len=85) :: linedata
END TYPE

end module mod_abundtypes

module mod_resultarrays

TYPE resultarray
        double precision :: NC_abund_CEL=0d0
        double precision :: cii_abund_CEL=0d0
        double precision :: ciii_abund_CEL=0d0
        double precision :: civ_abund_CEL=0d0
        double precision :: C_icf_CEL=0d0
        double precision :: C_abund_CEL=0d0
        double precision :: Nii_abund_CEL=0d0
        double precision :: Niii_abund_CEL=0d0
        double precision :: Niv_abund_CEL=0d0
        double precision :: Nv_abund_CEL=0d0
        double precision :: N_icf_CEL=0d0
        double precision :: N_abund_CEL=0d0
        double precision :: NO_abund_CEL=0d0
        double precision :: Oii_abund_CEL=0d0
        double precision :: Oiii_abund_CEL=0d0
        double precision :: Oiv_abund_CEL=0d0
        double precision :: O_icf_CEL=0d0
        double precision :: O_abund_CEL=0d0
        double precision :: Neii_abund_CEL=0d0
        double precision :: Neiii_abund_CEL=0d0
        double precision :: Neiv_abund_CEL=0d0
        double precision :: Nev_abund_CEL=0d0
        double precision :: Ne_icf_CEL=0d0
        double precision :: Ne_abund_CEL=0d0
        double precision :: Ariii_abund_CEL=0d0
        double precision :: Ariv_abund_CEL=0d0
        double precision :: Arv_abund_CEL=0d0
        double precision :: Ar_icf_CEL=0d0
        double precision :: Ar_abund_CEL=0d0
        double precision :: Sii_abund_CEL=0d0
        double precision :: Siii_abund_CEL=0d0
        double precision :: Cliii_abund_CEL=0d0
        double precision :: Cl_icf_CEL=0d0
        double precision :: Cl_abund_CEL=0d0
        double precision :: S_icf_CEL=0d0
        double precision :: S_abund_CEL=0d0
        double precision :: Hei_abund_ORL=0d0
        double precision :: Heii_abund_ORL=0d0
        double precision :: He_abund_ORL=0d0
        double precision :: Cii_abund_ORL=0d0
        double precision :: Ciii_abund_ORL=0d0
        double precision :: C_icf_ORL=0d0
        double precision :: C_abund_ORL=0d0
        double precision :: Nii_abund_ORL=0d0
        double precision :: Niii_abund_ORL=0d0
        double precision :: N_icf_ORL=0d0
        double precision :: N_abund_ORL=0d0
        double precision :: Oii_abund_ORL=0d0
        double precision :: O_icf_ORL=0d0
        double precision :: O_abund_ORL=0d0
        double precision :: Neii_abund_ORL=0d0
        double precision :: Ne_icf_ORL=0d0
        double precision :: Ne_abund_ORL=0d0
        double precision :: OII_density=0d0
        double precision :: SII_density=0d0
        double precision :: low_density=0d0
        double precision :: OII_temp=0d0
        double precision :: NII_temp=0d0
        double precision :: SII_temp=0d0
        double precision :: OI_temp=0d0
        double precision :: CI_temp=0d0
        double precision :: low_temp=0d0
        double precision :: ClIII_density=0d0
        double precision :: ArIV_density=0d0
        double precision :: CIII_density=0d0
        double precision :: OIII_IR_density=0d0
        double precision :: SIII_IR_density=0d0
        double precision :: ArIII_IR_density=0d0
        double precision :: NeIII_IR_density=0d0
        double precision :: med_density=0d0
        double precision :: OIII_temp=0d0
        double precision :: OIII_IR_temp=0d0
        double precision :: OIII_UV_temp=0d0
        double precision :: NeIII_temp=0d0
        double precision :: NeIII_IR_temp=0d0
        double precision :: ArIII_temp=0d0
        double precision :: SIII_temp=0d0
        double precision :: med_temp=0d0
        double precision :: NeIV_density=0d0
        double precision :: high_density=0d0
        double precision :: ArV_temp=0d0
        double precision :: NeV_temp=0d0
        double precision :: high_temp=0d0
        double precision :: mean_cHb=0d0
        double precision :: O_R23_upper=0d0
        double precision :: O_R23_lower=0d0
        double precision :: O_N2=0d0
        double precision :: O_O3N2=0d0
        double precision :: O_Ar3O3=0d0
        double precision :: O_S3O3=0d0
        double precision :: adf_O=0d0
        double precision :: adf_O2plus=0d0
        double precision :: adf_N=0d0
        double precision :: adf_N2plus=0d0
        double precision :: adf_C=0d0
        double precision :: adf_C2plus=0d0
        double precision :: adf_Ne=0d0
        double precision :: adf_Ne2plus=0d0
        double precision :: OII_density_ratio=0d0
        double precision :: SII_density_ratio=0d0
        double precision :: NII_temp_ratio=0d0
        double precision :: OII_temp_ratio=0d0
        double precision :: SII_temp_ratio=0d0
        double precision :: CI_temp_ratio=0d0
        double precision :: OI_temp_ratio=0d0
        double precision :: ClIII_density_ratio=0d0
        double precision :: ArIV_density_ratio=0d0
        double precision :: CIII_density_ratio=0d0
        double precision :: OIII_IR_density_ratio=0d0
        double precision :: ArIII_IR_density_ratio=0d0
        double precision :: SIII_IR_density_ratio=0d0
        double precision :: NeIII_IR_density_ratio=0d0
        double precision :: OIII_temp_ratio=0d0
        double precision :: NeIII_temp_ratio=0d0
        double precision :: ArIII_temp_ratio=0d0
        double precision :: SIII_temp_ratio=0d0
        double precision :: OIII_IR_temp_ratio=0d0
        double precision :: OIII_UV_temp_ratio=0d0
        double precision :: NeIII_IR_temp_ratio=0d0
        double precision :: NeIV_density_ratio=0d0
        double precision :: ArV_temp_ratio=0d0
        double precision :: NeV_temp_ratio=0d0
        double precision :: Bal_jump_temp=0d0
        double precision :: oii_v1_abund_orl=0d0
        double precision :: oii_v2_abund_orl=0d0
        double precision :: oii_v5_abund_orl=0d0
        double precision :: oii_v10_abund_orl=0d0
        double precision :: oii_v11_abund_orl=0d0
        double precision :: oii_v12_abund_orl=0d0
        double precision :: oii_v19_abund_orl=0d0
        double precision :: oii_v20_abund_orl=0d0
        double precision :: oii_v25_abund_orl=0d0
        double precision :: oii_v28_abund_orl=0d0
        double precision :: oii_v33_abund_orl=0d0
        double precision :: oii_3d4f_abund_orl=0d0
        double precision :: nii_v3_abund_orl=0d0
        double precision :: nii_v5_abund_orl=0d0
        double precision :: nii_v8_abund_orl=0d0
        double precision :: nii_v12_abund_orl =0d0
        double precision :: nii_v20_abund_orl =0d0
        double precision :: nii_v28_abund_orl=0d0
        double precision :: nii_3d4f_abund_orl=0d0
end type

end module mod_resultarrays

module mod_atomicdata

type atomic_data

        integer :: NTEMPS
        integer :: NLEVS
        integer :: irats
        character(len=10) :: ion
        real(kind=8),allocatable :: temps(:)
        real(kind=8),allocatable :: roott(:)
        character(len=20),allocatable :: labels(:)
        real(kind=8),allocatable :: G(:)
        real(kind=8),allocatable :: waveno(:)
        real(kind=8),allocatable :: A_coeffs(:,:)
        real(kind=8),allocatable :: col_str(:,:,:)

end type

end module mod_atomicdata
