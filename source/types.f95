module mod_abundtypes

implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)

TYPE line
        CHARACTER(len=11) :: name
        CHARACTER(len=20) :: ion
        real(kind=dp) :: wavelength
        real(kind=dp) :: wavelength_observed !for use with output from ALFA
        real(kind=dp) :: intensity
        real(kind=dp) :: int_err
        real(kind=dp) :: int_dered
        real(kind=dp) :: blend_intensity !intensity for blended lines copied to this, set to zero for abundance calculations, then copied back for linelist
        real(kind=dp) :: blend_int_dered
        real(kind=dp) :: blend_int_err
        real(kind=dp) :: freq
        CHARACTER(len=20) :: transition
        real(kind=dp) :: abundance
        CHARACTER(len=4) :: zone !high, medium or low ionisation zone
        INTEGER :: location !for CELs and He, this indicates the position of the line in the main linelist array, so that when its abundance is calculated it can be copied back into the linelist array for easy outputting
        CHARACTER(len=15) :: latextext ! ion name in latex format
        CHARACTER(len=85) :: linedata
END TYPE

end module mod_abundtypes

module mod_resultarrays

implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)

TYPE resultarray
        real(kind=dp) :: NC_abund_CEL=0d0
        real(kind=dp) :: cii_abund_CEL=0d0
        real(kind=dp) :: ciii_abund_CEL=0d0
        real(kind=dp) :: civ_abund_CEL=0d0
        real(kind=dp) :: C_icf_CEL=0d0
        real(kind=dp) :: C_abund_CEL=0d0
        real(kind=dp) :: Nii_abund_CEL=0d0
        real(kind=dp) :: Niii_abund_CEL=0d0
        real(kind=dp) :: Niv_abund_CEL=0d0
        real(kind=dp) :: Nv_abund_CEL=0d0
        real(kind=dp) :: N_icf_CEL=0d0
        real(kind=dp) :: N_abund_CEL=0d0
        real(kind=dp) :: NO_abund_CEL=0d0
        real(kind=dp) :: Oii_abund_CEL=0d0
        real(kind=dp) :: Oiii_abund_CEL=0d0
        real(kind=dp) :: Oiv_abund_CEL=0d0
        real(kind=dp) :: O_icf_CEL=0d0
        real(kind=dp) :: O_abund_CEL=0d0
        real(kind=dp) :: Neii_abund_CEL=0d0
        real(kind=dp) :: Neiii_abund_CEL=0d0
        real(kind=dp) :: Neiv_abund_CEL=0d0
        real(kind=dp) :: Nev_abund_CEL=0d0
        real(kind=dp) :: Ne_icf_CEL=0d0
        real(kind=dp) :: Ne_abund_CEL=0d0
        real(kind=dp) :: Ariii_abund_CEL=0d0
        real(kind=dp) :: Ariv_abund_CEL=0d0
        real(kind=dp) :: Arv_abund_CEL=0d0
        real(kind=dp) :: Ar_icf_CEL=0d0
        real(kind=dp) :: Ar_abund_CEL=0d0
        real(kind=dp) :: Sii_abund_CEL=0d0
        real(kind=dp) :: Siii_abund_CEL=0d0
        real(kind=dp) :: Cliii_abund_CEL=0d0
        real(kind=dp) :: Cl_icf_CEL=0d0
        real(kind=dp) :: Cl_abund_CEL=0d0
        real(kind=dp) :: S_icf_CEL=0d0
        real(kind=dp) :: S_abund_CEL=0d0
        real(kind=dp) :: Hei_abund_ORL=0d0
        real(kind=dp) :: Heii_abund_ORL=0d0
        real(kind=dp) :: He_abund_ORL=0d0
        real(kind=dp) :: Cii_abund_ORL=0d0
        real(kind=dp) :: Ciii_abund_ORL=0d0
        real(kind=dp) :: C_icf_ORL=0d0
        real(kind=dp) :: C_abund_ORL=0d0
        real(kind=dp) :: Nii_abund_ORL=0d0
        real(kind=dp) :: Niii_abund_ORL=0d0
        real(kind=dp) :: N_icf_ORL=0d0
        real(kind=dp) :: N_abund_ORL=0d0
        real(kind=dp) :: Oii_abund_ORL=0d0
        real(kind=dp) :: O_icf_ORL=0d0
        real(kind=dp) :: O_abund_ORL=0d0
        real(kind=dp) :: Neii_abund_ORL=0d0
        real(kind=dp) :: Ne_icf_ORL=0d0
        real(kind=dp) :: Ne_abund_ORL=0d0
        real(kind=dp) :: OII_density=0d0
        real(kind=dp) :: SII_density=0d0
        real(kind=dp) :: low_density=0d0
        real(kind=dp) :: OII_temp=0d0
        real(kind=dp) :: NII_temp=0d0
        real(kind=dp) :: SII_temp=0d0
        real(kind=dp) :: OI_temp=0d0
        real(kind=dp) :: CI_temp=0d0
        real(kind=dp) :: low_temp=0d0
        real(kind=dp) :: ClIII_density=0d0
        real(kind=dp) :: ArIV_density=0d0
        real(kind=dp) :: CIII_density=0d0
        real(kind=dp) :: OIII_IR_density=0d0
        real(kind=dp) :: SIII_IR_density=0d0
        real(kind=dp) :: ArIII_IR_density=0d0
        real(kind=dp) :: NeIII_IR_density=0d0
        real(kind=dp) :: med_density=0d0
        real(kind=dp) :: OIII_temp=0d0
        real(kind=dp) :: OIII_IR_temp=0d0
        real(kind=dp) :: OIII_UV_temp=0d0
        real(kind=dp) :: NeIII_temp=0d0
        real(kind=dp) :: NeIII_IR_temp=0d0
        real(kind=dp) :: ArIII_temp=0d0
        real(kind=dp) :: SIII_temp=0d0
        real(kind=dp) :: med_temp=0d0
        real(kind=dp) :: NeIV_density=0d0
        real(kind=dp) :: high_density=0d0
        real(kind=dp) :: ArV_temp=0d0
        real(kind=dp) :: NeV_temp=0d0
        real(kind=dp) :: high_temp=0d0
        real(kind=dp) :: mean_cHb=0d0
        real(kind=dp) :: O_R23_upper=0d0
        real(kind=dp) :: O_R23_lower=0d0
        real(kind=dp) :: O_N2=0d0
        real(kind=dp) :: O_O3N2=0d0
        real(kind=dp) :: O_Ar3O3=0d0
        real(kind=dp) :: O_S3O3=0d0
        real(kind=dp) :: adf_O=0d0
        real(kind=dp) :: adf_O2plus=0d0
        real(kind=dp) :: adf_N=0d0
        real(kind=dp) :: adf_N2plus=0d0
        real(kind=dp) :: adf_C=0d0
        real(kind=dp) :: adf_C2plus=0d0
        real(kind=dp) :: adf_Ne=0d0
        real(kind=dp) :: adf_Ne2plus=0d0
        real(kind=dp) :: OII_density_ratio=0d0
        real(kind=dp) :: SII_density_ratio=0d0
        real(kind=dp) :: NII_temp_ratio=0d0
        real(kind=dp) :: OII_temp_ratio=0d0
        real(kind=dp) :: SII_temp_ratio=0d0
        real(kind=dp) :: CI_temp_ratio=0d0
        real(kind=dp) :: OI_temp_ratio=0d0
        real(kind=dp) :: ClIII_density_ratio=0d0
        real(kind=dp) :: ArIV_density_ratio=0d0
        real(kind=dp) :: CIII_density_ratio=0d0
        real(kind=dp) :: OIII_IR_density_ratio=0d0
        real(kind=dp) :: ArIII_IR_density_ratio=0d0
        real(kind=dp) :: SIII_IR_density_ratio=0d0
        real(kind=dp) :: NeIII_IR_density_ratio=0d0
        real(kind=dp) :: OIII_temp_ratio=0d0
        real(kind=dp) :: NeIII_temp_ratio=0d0
        real(kind=dp) :: ArIII_temp_ratio=0d0
        real(kind=dp) :: SIII_temp_ratio=0d0
        real(kind=dp) :: OIII_IR_temp_ratio=0d0
        real(kind=dp) :: OIII_UV_temp_ratio=0d0
        real(kind=dp) :: NeIII_IR_temp_ratio=0d0
        real(kind=dp) :: NeIV_density_ratio=0d0
        real(kind=dp) :: ArV_temp_ratio=0d0
        real(kind=dp) :: NeV_temp_ratio=0d0
        real(kind=dp) :: Bal_jump_temp=0d0
        real(kind=dp) :: oii_v1_abund_orl=0d0
        real(kind=dp) :: oii_v2_abund_orl=0d0
        real(kind=dp) :: oii_v5_abund_orl=0d0
        real(kind=dp) :: oii_v10_abund_orl=0d0
        real(kind=dp) :: oii_v11_abund_orl=0d0
        real(kind=dp) :: oii_v12_abund_orl=0d0
        real(kind=dp) :: oii_v19_abund_orl=0d0
        real(kind=dp) :: oii_v20_abund_orl=0d0
        real(kind=dp) :: oii_v25_abund_orl=0d0
        real(kind=dp) :: oii_v28_abund_orl=0d0
        real(kind=dp) :: oii_v33_abund_orl=0d0
        real(kind=dp) :: oii_3d4f_abund_orl=0d0
        real(kind=dp) :: nii_v3_abund_orl=0d0
        real(kind=dp) :: nii_v5_abund_orl=0d0
        real(kind=dp) :: nii_v8_abund_orl=0d0
        real(kind=dp) :: nii_v12_abund_orl =0d0
        real(kind=dp) :: nii_v20_abund_orl =0d0
        real(kind=dp) :: nii_v28_abund_orl=0d0
        real(kind=dp) :: nii_3d4f_abund_orl=0d0
        real(kind=dp) :: oii_te=0.d0
        real(kind=dp) :: oii_ne=0.d0
        real(kind=dp) :: ratio_5876_4471=0.d0
        real(kind=dp) :: ratio_6678_4471=0.d0
        real(kind=dp) :: Te_5876_4471=0.d0
        real(kind=dp) :: Te_6678_4471=0.d0
        real(kind=dp) :: nii5754recCEL=0.d0
        real(kind=dp) :: oii7325recCEL=0.d0
        real(kind=dp) :: oiii4363recCEL=0.d0
        real(kind=dp) :: nii5754recRL=0.d0
        real(kind=dp) :: oii7325recRL=0.d0
        real(kind=dp) :: oiii4363recRL=0.d0
        real(kind=dp) :: balmerdec_density=0.d0

end type

end module mod_resultarrays

module mod_atomicdata

implicit none
private :: dp
integer, parameter :: dp = kind(1.d0)

type atomic_data

        integer :: NTEMPS
        integer :: NLEVS
        integer :: irats
        character(len=10) :: ion
        real(kind=dp),allocatable :: temps(:)
        real(kind=dp),allocatable :: roott(:)
        character(len=20),allocatable :: labels(:)
        integer,allocatable :: G(:)
        real(kind=dp),allocatable :: waveno(:)
        real(kind=dp),allocatable :: A_coeffs(:,:)
        real(kind=dp),allocatable :: col_str(:,:,:)

end type

end module mod_atomicdata
