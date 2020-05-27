!rec_lines.f90, routines to calculate ionic abundances from recombination lines
!(C) Roger Wesson, Dave Stock, Peter Scicluna
! based originally on a conversion of MIDAS script Roii.prg, written by XWL, to F90
! RW May 2009

module mod_recombination_lines
use mod_globals

      implicit none

      type oiiRL
            character(len=1) :: Hyb
            character(len=1) :: n_E1
            character(len=1) :: n_E1GA
            character(len=1) :: n_E2
            character(len=1) :: n_E2GA
            character(len=1) :: n_g1
            character(len=1) :: n_g2
            character(len=1) :: Rem1
            character(len=1) :: Rem2
            character(len=1) :: Rem3
            character(len=1) :: Rem4
            character(len=3) :: q_gf1
            character(len=3) :: q_gf2
            character(len=7) :: Mult
            character(len=9) :: Term1
            character(len=9) :: Term2
            integer :: g1
            integer :: g2
            integer :: ION
            real(kind=dp) :: Wave
            real(kind=dp) :: E1
            real(kind=dp) :: E2
            real(kind=dp) :: Em
            real(kind=dp) :: Int
            real(kind=dp) :: Br_A
            real(kind=dp) :: Br_B
            real(kind=dp) :: Br_C
            real(kind=dp) :: gf1
            real(kind=dp) :: gf2
            real(kind=dp) :: Obs
            real(kind=dp) :: abundance
      end type oiiRL

      type(oiiRL), dimension(:), allocatable :: oiiRLs

      type oiiRL_s2017
        real(kind=dp) :: Wave
        character(len=8) :: Mult
        character(len=9) :: Term1
        character(len=9) :: Term2
        real(kind=dp) :: Int
        real(kind=dp) :: Obs
        real(kind=dp) :: abundance
        real(kind=dp), dimension(25,16) :: ems
        real(kind=dp) :: em_interpolated
      end type oiiRL_s2017

      type(oiiRL_s2017), dimension(276) :: oiiRLs_s2017

      real(kind=dp), dimension(36,9), target :: oii_coefficients
      real(kind=dp), dimension(:), pointer :: oii_A_4f, oii_A_3d4F, oii_A_3d4D, oii_B_3d4D, oii_C_3d4D, oii_A_3d2F, oii_B_3d2F, oii_C_3d2F, oii_A_3d2D, oii_C_3d2D, oii_A_3d2P, oii_C_3d2P, oii_A_3p4D, oii_B_3p4D, oii_A_3p4P, oii_B_3p4P, oii_A_3p4S, oii_B_3p4S, oii_A_3p2D, oii_C_3p2D, oii_A_3p2P, oii_C_3p2P, oii_A_3p2S, oii_C_3p2S, oii_A_4f_low, oii_A_3d4F_low, oii_B_3d4D_low, oii_A_3d2F_low, oii_A_3d2D_low, oii_A_3d2P_low, oii_A_3p4D_low, oii_A_3p4P_low, oii_A_3p4S_low, oii_A_3p2D_low, oii_A_3p2P_low, oii_A_3p2S_low

      real(kind=dp), dimension(27,6), target :: nii_coefficients
      real(kind=dp), dimension(:), pointer :: nii_V3, nii_V4, nii_V5, nii_V8, nii_V12, nii_V13, nii_V15, nii_V17, nii_V19, nii_V20, nii_V21, nii_V22, nii_V24, nii_V26, nii_V28, nii_V29, nii_V30, nii_V31, nii_V36, nii_V39, nii_V58, nii_V48, nii_V55, nii_V43, nii_V61, nii_V39b, nii_V58a

      type niiRL
            character(len=1) :: Hyb
            character(len=1) :: n_E1
            character(len=1) :: n_E1GA
            character(len=1) :: n_E2
            character(len=1) :: n_E2GA
            character(len=1) :: n_g1
            character(len=1) :: n_g2
            character(len=1) :: Rem1
            character(len=1) :: Rem2
            character(len=1) :: Rem3
            character(len=1) :: Rem4
            character(len=3) :: q_gf1
            character(len=3) :: q_gf2
            character(len=7) :: Mult
            character(len=9) :: Term1
            character(len=9) :: Term2
            integer :: g1
            integer :: g2
            integer :: ION
            real(kind=dp) :: Wave
            real(kind=dp) :: E1
            real(kind=dp) :: E2
            real(kind=dp) :: Em
            real(kind=dp) :: Int
            real(kind=dp) :: Br_LS
            real(kind=dp) :: gf1
            real(kind=dp) :: gf2
            real(kind=dp) :: Obs
            real(kind=dp) :: abundance
      end type niiRL

      type(niiRL), dimension(:),allocatable :: niiRLs

      type ciiRL
            real(kind=dp) :: Wave
            real(kind=dp) :: a
            real(kind=dp) :: b
            real(kind=dp) :: c
            real(kind=dp) :: d
            real(kind=dp) :: f
            real(kind=dp) :: aeff
            real(kind=dp) :: Int
            real(kind=dp) :: Obs
            real(kind=dp) :: abundance
      end type ciiRL

      type(ciiRL), dimension(:),allocatable :: ciiRLs

      type neiiRL
            real(kind=dp) :: Wave
            real(kind=dp) :: a
            real(kind=dp) :: b
            real(kind=dp) :: c
            real(kind=dp) :: d
            real(kind=dp) :: f
            real(kind=dp) :: Br
            real(kind=dp) :: aeff
            real(kind=dp) :: Int
            real(kind=dp) :: Obs
            real(kind=dp) :: abundance
      end type neiiRL

      type(neiiRL), dimension(:),allocatable :: neiiRLs

      type xiiiRL
            character(len=3) :: Ion
            real(kind=dp) :: Wave
            real(kind=dp) :: a
            real(kind=dp) :: b
            real(kind=dp) :: c
            real(kind=dp) :: d
            real(kind=dp) :: Br
            real(kind=dp) :: aeff
            real(kind=dp) :: Int
            real(kind=dp) :: Obs
            real(kind=dp) :: abundance
      end type xiiiRL

      type(xiiiRL), dimension(:),allocatable :: xiiiRLs

      contains

subroutine read_orl_data

        implicit none
        integer :: i, nlines

        ! read in OII data

!debugging
#ifdef CO
        print *,"subroutine: read_orl_data"
#endif

            301 format (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, F7.4,     &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
            open(201, file=trim(PREFIX)//"/share/neat/Roii.dat", status='old')
            read(201,*) nlines
            allocate(oiiRLs(nlines))
            oiiRLs%Int = 0.d0
            oiiRLs%Obs=0.d0
            oiiRLs%abundance=0.d0
            do i = 1,nlines
            read(201,301) oiiRLs(i)%ion, oiiRLs(i)%Wave, oiiRLs(i)%Hyb, &
     &oiiRLs(i)%Rem1, oiiRLs(i)%Rem2, oiiRLs(i)%Rem3, oiiRLs(i)%Rem4,   &
     &oiiRLs(i)%gf1, oiiRLs(i)%q_gf1, oiiRLs(i)%gf2, oiiRLs(i)%q_gf2,   &
     &oiiRLs(i)%Mult, oiiRLs(i)%E1, oiiRLs(i)%n_E1, oiiRLs(i)%n_E1GA,   &
     &oiiRLs(i)%g1, oiiRLs(i)%n_g1, oiiRLs(i)%Term1, oiiRLs(i)%E2,      &
     &oiiRLs(i)%n_E2, oiiRLs(i)%n_E2GA, oiiRLs(i)%g2, oiiRLs(i)%n_g2,   &
     &oiiRLs(i)%Term2, oiiRLs(i)%Br_A, oiiRLs(i)%Br_B, oiiRLs(i)%Br_C
            enddo
      close(201)

     ! coefficients from LSBC 1995, S94
     ! array consists of coefficients a2, a4, a5, a6, b, c and d, plus calculated values a and aeff
     ! an is the coefficient a at log(ne)=n. this will be interpolated and stored in variable a.
     ! the interpolated a is then used with the other coefficients to calculate aeff
     ! replace with reading from file at some point

      oii_coefficients(1,:) = (/0.236d0,0.232d0,0.228d0,0.222d0,-0.92009d0,0.15526d0,0.03442d0,0.d0,0.d0/)
      oii_coefficients(2,:) = (/0.876d0,0.876d0,0.877d0,0.880d0,-0.73465d0,0.13689d0,0.06220d0,0.d0,0.d0/)
!      oii_coefficients(3,:) = (/0.727d0,0.726d0,0.725d0,0.726d0,-0.73465d0,0.13689d0,0.06220d0,0.d0,0.d0/)
      oii_coefficients(4,:) = (/0.747d0,0.745d0,0.744d0,0.745d0,-0.74621d0,0.15710d0,0.07059d0,0.d0,0.d0/)
!      oii_coefficients(5,:) = (/0.769d0,0.767d0,0.766d0,0.766d0,-0.74621d0,0.15710d0,0.07059d0,0.d0,0.d0/)
      oii_coefficients(6,:) = (/0.727d0,0.726d0,0.725d0,0.726d0,-0.74621d0,0.15710d0,0.07059d0,0.d0,0.d0/)
!      oii_coefficients(7,:) = (/0.747d0,0.745d0,0.744d0,0.745d0,-0.74621d0,0.15710d0,0.07059d0,0.d0,0.d0/)
!      oii_coefficients(8,:) = (/0.769d0,0.767d0,0.766d0,0.766d0,-0.74621d0,0.15710d0,0.07059d0,0.d0,0.d0/)
      oii_coefficients(9,:) = (/0.603d0,0.601d0,0.600d0,0.599d0,-0.79533d0,0.15314d0,0.05322d0,0.d0,0.d0/)
!      oii_coefficients(10,:) = (/0.620d0,0.618d0,0.616d0,0.615d0,-0.79533d0,0.15314d0,0.05322d0,0.d0,0.d0/)
      oii_coefficients(11,:) = (/0.526d0,0.524d0,0.523d0,0.524d0,-0.78448d0,0.13681d0,0.05608d0,0.d0,0.d0/)
!      oii_coefficients(1,:) = (/0.538d0,0.536d0,0.535d0,0.536d0,-0.78448d0,0.13681d0,0.05608d0,0.d0,0.d0/)
      oii_coefficients(13,:) = (/34.7d0,34.9d0,35.1d0,35.0d0,-0.749d0,0.023d0,0.074d0,0.d0,0.d0/)
      oii_coefficients(14,:) = (/36.0d0,36.2d0,36.4d0,36.3d0,-0.736d0,0.033d0,0.077d0,0.d0,0.d0/)
      oii_coefficients(15,:) = (/10.4d0,10.4d0,10.5d0,10.4d0,-0.721d0,0.073d0,0.072d0,0.d0,0.d0/)
      oii_coefficients(16,:) = (/14.6d0,14.6d0,14.7d0,14.6d0,-0.732d0,0.081d0,0.066d0,0.d0,0.d0/)
      oii_coefficients(17,:) = (/0.90d0,0.90d0,0.90d0,1.00d0,-0.485d0,-0.047d0,0.140d0,0.d0,0.d0/)
      oii_coefficients(18,:) = (/4.80d0,4.90d0,4.90d0,4.90d0,-0.730d0,-0.003d0,0.057d0,0.d0,0.d0/)
      oii_coefficients(19,:) = (/2.40d0,2.40d0,2.50d0,2.60d0,-0.550d0,-0.051d0,0.178d0,0.d0,0.d0/)
      oii_coefficients(20,:) = (/14.5d0,14.6d0,14.5d0,14.3d0,-0.736d0,0.068d0,0.066d0,0.d0,0.d0/)
      oii_coefficients(21,:) = (/1.10d0,1.20d0,1.20d0,1.20d0,-0.523d0,-0.044d0,0.173d0,0.d0,0.d0/)
      oii_coefficients(22,:) = (/1.30d0,1.40d0,1.40d0,1.40d0,-0.565d0,-0.042d0,0.158d0,0.d0,0.d0/)
      oii_coefficients(23,:) = (/0.40d0,0.40d0,0.40d0,0.40d0,-0.461d0,-0.083d0,0.287d0,0.d0,0.d0/)
      oii_coefficients(24,:) = (/0.50d0,0.50d0,0.50d0,0.60d0,-0.547d0,-0.074d0,0.244d0,0.d0,0.d0/)

      oii_coefficients(25,:) = (/0.236d0,0.236d0,0.236d0,0.236d0,-1.07552d0,-0.04843d0,0.d0,0.d0,0.d0/)
      oii_coefficients(26,:) = (/0.878d0,0.878d0,0.878d0,0.878d0,-0.86175d0,-0.02470d0,0.d0,0.d0,0.d0/)
      oii_coefficients(27,:) = (/0.747d0,0.747d0,0.747d0,0.747d0,-0.89382d0,-0.02906d0,0.d0,0.d0,0.d0/)
      oii_coefficients(28,:) = (/0.747d0,0.747d0,0.747d0,0.747d0,-0.89382d0,-0.02906d0,0.d0,0.d0,0.d0/)
      oii_coefficients(29,:) = (/0.603d0,0.603d0,0.603d0,0.603d0,-0.94025d0,-0.03467d0,0.d0,0.d0,0.d0/)
      oii_coefficients(30,:) = (/0.526d0,0.526d0,0.526d0,0.526d0,-0.91758d0,-0.03120d0,0.d0,0.d0,0.d0/)
      oii_coefficients(31,:) = (/36.288d0,36.288d0,36.288d0,36.288d0,-0.75421d0,0.02883d0,0.01213d0,0.d0,0.d0/)
      oii_coefficients(32,:) = (/14.656d0,14.656d0,14.656d0,14.656d0,-0.80449d0,0.00018d0,0.00517d0,0.d0,0.d0/)
      oii_coefficients(33,:) = (/4.8340d0,4.8340d0,4.8340d0,4.8340d0,-0.71947d0,0.02544d0,0.00936d0,0.d0,0.d0/)
      oii_coefficients(34,:) = (/2.3616d0,2.3616d0,2.3616d0,2.3616d0,-0.46263d0,0.14697d0,0.03856d0,0.d0,0.d0/)
      oii_coefficients(35,:) = (/1.1198d0,1.1198d0,1.1198d0,1.1198d0,-0.44147d0,0.13837d0,0.03191d0,0.d0,0.d0/)
      oii_coefficients(36,:) = (/0.3922d0,0.3922d0,0.3922d0,0.3922d0,-0.35043d0,0.26366d0,0.06666d0,0.d0,0.d0/)

!define pointers so that we can refer to the coefficients by state as well as processing all coefficients at once

      oii_A_4f => oii_coefficients(1,:)
    oii_A_3d4F => oii_coefficients(2,:)
    oii_A_3d4D => oii_coefficients(3,:)
    oii_B_3d4D => oii_coefficients(4,:)
    oii_C_3d4D => oii_coefficients(5,:)
    oii_A_3d2F => oii_coefficients(6,:)
    oii_B_3d2F => oii_coefficients(7,:)
    oii_C_3d2F => oii_coefficients(8,:)
    oii_A_3d2D => oii_coefficients(9,:)
    oii_C_3d2D => oii_coefficients(10,:)
    oii_A_3d2P => oii_coefficients(11,:)
    oii_C_3d2P => oii_coefficients(12,:)
    oii_A_3p4D => oii_coefficients(13,:)
    oii_B_3p4D => oii_coefficients(14,:)
    oii_A_3p4P => oii_coefficients(15,:)
    oii_B_3p4P => oii_coefficients(16,:)
    oii_A_3p4S => oii_coefficients(17,:)
    oii_B_3p4S => oii_coefficients(18,:)
    oii_A_3p2D => oii_coefficients(19,:)
    oii_C_3p2D => oii_coefficients(20,:)
    oii_A_3p2P => oii_coefficients(21,:)
    oii_C_3p2P => oii_coefficients(22,:)
    oii_A_3p2S => oii_coefficients(23,:)
    oii_C_3p2S => oii_coefficients(24,:)
  oii_A_4f_low => oii_coefficients(25,:)
oii_A_3d4F_low => oii_coefficients(26,:)
oii_B_3d4D_low => oii_coefficients(27,:)
oii_A_3d2F_low => oii_coefficients(28,:)
oii_A_3d2D_low => oii_coefficients(29,:)
oii_A_3d2P_low => oii_coefficients(30,:)
oii_A_3p4D_low => oii_coefficients(31,:)
oii_A_3p4P_low => oii_coefficients(32,:)
oii_A_3p4S_low => oii_coefficients(33,:)
oii_A_3p2D_low => oii_coefficients(34,:)
oii_A_3p2P_low => oii_coefficients(35,:)
oii_A_3p2S_low => oii_coefficients(36,:)

     ! read in NII data

            302 format (I5, 1X, F9.4, 1X, A1, A1, A1, A1, A1, 1X, F7.4, &
     & 1X, A3, 1X, F7.4, 1X, A3, 1X, A7, 3X, F11.4, A1, A1, 1X, I2, &
     &1X, A1, 1X, A9, 1X, F13.4, 1X, A1, A1, 1X, I2, 1X, A1, 1X, A9, 1X,&
     & F7.4, 1X, F7.4, 1X, F7.4)!, 1X, E10.4, 1X, E10.4, 1X)
            open(201, file=trim(PREFIX)//"/share/neat/Rnii.dat", status='old')
            read(201,*) nlines
            allocate(niiRLs(nlines))
            niiRLs%Int = 0.d0
            niiRLs%Obs=0.d0
            niiRLs%abundance=0.d0
            do i = 1,nlines
            read(201,302) niiRLs(i)%ION, niiRLs(i)%Wave, niiRLs(i)%Hyb, &
     &niiRLs(i)%Rem1, niiRLs(i)%Rem2, niiRLs(i)%Rem3, niiRLs(i)%Rem4,   &
     &niiRLs(i)%gf1, niiRLs(i)%q_gf1, niiRLs(i)%gf2, niiRLs(i)%q_gf2,   &
     &niiRLs(i)%Mult, niiRLs(i)%E1, niiRLs(i)%n_E1, niiRLs(i)%n_E1GA,   &
     &niiRLs(i)%g1, niiRLs(i)%n_g1, niiRLs(i)%Term1, niiRLs(i)%E2,      &
     &niiRLs(i)%n_E2, niiRLs(i)%n_E2GA, niiRLs(i)%g2, niiRLs(i)%n_g2,   &
     &niiRLs(i)%Term2, niiRLs(i)%Br_LS
            enddo
      close(201)

! arrays of NII fitting coefficients, consisting of a, b, c, d, Br, and aeff

      nii_coefficients(1,:) = (/-12.7289d0, -0.689816d0, 0.022005d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(2,:) = (/-13.8161d0, -0.778606d0, -0.028944d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(3,:) = (/-13.0765d0, -0.734594d0, -0.0251909d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(4,:) = (/-14.1211d0, -0.608107d0, 0.0362301d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(5,:) = (/-13.7473d0, -0.509595d0, 0.0255685d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(6,:) = (/-14.3753d0, -0.515547d0, 0.0100966d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(7,:) = (/-14.3932d0, -0.887946d0, -0.0525855d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(8,:) = (/-15.0052d0, -0.89811d0, -0.0581789d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(9,:) = (/-12.6183d0, -0.840727d0, -0.0229685d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(10,:) = (/-13.3184d0, -0.884034d0, -0.0512093d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(11,:) = (/-14.5113d0, -0.87792d0, -0.0552785d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(12,:) = (/-14.1305d0, -0.487037d0, 0.0354135d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(13,:) = (/-13.3527d0, -0.878224d0, -0.0557112d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(14,:) = (/-14.9628d0, -0.486746d0, 0.0358261d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(15,:) = (/-13.0871d0, -0.883624d0, -0.0506882d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(16,:) = (/-13.5581d0, -0.878488d0, -0.0557583d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(17,:) = (/-14.3521d0, -0.487527d0, 0.0355516d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(18,:) = (/-15.0026d0, -0.923093d0, -0.0588371d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(19,:) = (/-13.8636d0, -0.569144d0, 0.0068655d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(20,:) = (/-13.035d0, -1.12035d0, -0.10642d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(21,:) = (/-13.5484d0, -1.11909d0, -0.105123d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(22,:) = (/-13.2548d0, -1.12902d0, -0.110368d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(23,:) = (/-13.5656d0, -1.11989d0, -0.105818d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(24,:) = (/-13.7426d0, -1.13351d0, -0.111146d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(25,:) = (/-13.7373d0, -1.12695d0, -0.108158d0, 0.d0, 0.d0, 0.d0 /)
      nii_coefficients(26,:) = (/0.108d0, -0.754d0, 2.587d0, 0.719d0, 0.350d0, 0.d0 /)
      nii_coefficients(27,:) = (/0.326d0, -0.754d0, 2.587d0, 0.719d0, 0.074d0, 0.d0 /)

        nii_V3 => nii_coefficients(1,:)
        nii_V4 => nii_coefficients(2,:)
        nii_V5 => nii_coefficients(3,:)
        nii_V8 => nii_coefficients(4,:)
       nii_V12 => nii_coefficients(5,:)
       nii_V13 => nii_coefficients(6,:)
       nii_V15 => nii_coefficients(7,:)
       nii_V17 => nii_coefficients(8,:)
       nii_V19 => nii_coefficients(9,:)
       nii_V20 => nii_coefficients(10,:)
       nii_V21 => nii_coefficients(11,:)
       nii_V22 => nii_coefficients(12,:)
       nii_V24 => nii_coefficients(13,:)
       nii_V26 => nii_coefficients(14,:)
       nii_V28 => nii_coefficients(15,:)
       nii_V29 => nii_coefficients(16,:)
       nii_V30 => nii_coefficients(17,:)
       nii_V31 => nii_coefficients(18,:)
       nii_V36 => nii_coefficients(19,:)
       nii_V39 => nii_coefficients(20,:)
       nii_V58 => nii_coefficients(21,:)
       nii_V48 => nii_coefficients(22,:)
       nii_V55 => nii_coefficients(23,:)
       nii_V43 => nii_coefficients(24,:)
       nii_V61 => nii_coefficients(25,:)
      nii_V39b => nii_coefficients(26,:)
      nii_V58a => nii_coefficients(27,:)

!define the pointers

! read in CII data

      303 format (F7.2, 1X, F6.4, 1X, F7.4, 1X, F7.4, 1X, F7.4, 1X, F7.4)
      open(201, file=trim(PREFIX)//"/share/neat/Rcii.dat", status='old')
      read(201,*) nlines
      allocate(ciiRLs(nlines))
      ciiRLs%Int = 0.d0
      ciiRLs%Obs=0.d0
      ciiRLs%abundance=0.d0
      do i = 1,nlines
        read(201,303) ciiRLs(i)%Wave, ciiRLs(i)%a, ciiRLs(i)%b, &
        & ciiRLs(i)%c, ciiRLs(i)%d, ciiRLs(i)%f
      enddo
      close(201)

          ! read in NeII data

      304 format (F7.2, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F6.3, 1X, F7.4, 1X, F6.3)
      open(201, file=trim(PREFIX)//"/share/neat/Rneii.dat", status='old')
      read(201,*) nlines
      allocate(neiiRLs(nlines))
      neiiRLs%Int = 0.d0
      neiiRLs%Obs=0.d0
      neiiRLs%abundance=0.d0
      do i = 1,nlines
        read(201,304) neiiRLs(i)%Wave, neiiRLs(i)%a, neiiRLs(i)%b, &
        & neiiRLs(i)%c, neiiRLs(i)%d, neiiRLs(i)%f, neiiRLs(i)%Br
      enddo
      close(201)

        ! read in XIII data

      305 format (A3,1X,F7.2, 1X, F5.3, 1X, F6.3, 1X, F5.3, 1X, F5.3, 1X, F5.4)
      open(201, file=trim(PREFIX)//"/share/neat/Rxiii.dat", status='old')
      read(201,*) nlines
      allocate(xiiiRLs(nlines))
      xiiiRLs%Int = 0.d0
      xiiiRLs%Obs=0.d0
      xiiiRLs%abundance=0.d0
      do i = 1,nlines
        read(201,305) xiiiRLs(i)%ion, xiiiRLs(i)%Wave, xiiiRLs(i)%a, &
        & xiiiRLs(i)%b, xiiiRLs(i)%c, xiiiRLs(i)%d, xiiiRLs(i)%Br
      enddo
      close(201)

end subroutine read_orl_data

subroutine oii_rec_lines(te,ne,abund,oiiRLs,aeff_hb,em_hb)

      implicit none
      real(kind=dp) :: aeff, aeff_hb, em_hb, Te, Ne, abund, tered

      type(oiiRL), dimension(:) :: oiiRLs

!debugging
#ifdef CO
        print *,"subroutine: oii_rec_lines"
#endif

! interpolate the a values

      if (log10(ne) .le. 2) then
        oii_coefficients(:,8)= oii_coefficients(:,1)
      elseif (log10(ne) .gt. 2 .and. log10(ne) .le. 4) then
        oii_coefficients(:,8)= oii_coefficients(:,1) + (oii_coefficients(:,2) - oii_coefficients(:,1)) / 2. * (log10(ne) - 2.)
      elseif (log10(ne) .gt. 4 .and. log10(ne) .le. 5) then
        oii_coefficients(:,8)= oii_coefficients(:,2) + (oii_coefficients(:,3) - oii_coefficients(:,2)) * (log10(ne) - 2.)
      elseif (log10(ne) .gt. 6 .and. log10(ne) .le. 6) then
        oii_coefficients(:,8)= oii_coefficients(:,3) + (oii_coefficients(:,4) - oii_coefficients(:,3)) * (log10(ne) - 2.)
      else
        oii_coefficients(:,8)= oii_coefficients(:,4)
      endif

! calculate the aeffs.  coefficients

      tered=te/10000.
      oii_coefficients(1:24,9)=1.e-14 * (oii_coefficients(1:24,8) * tered**oii_coefficients(1:24,5) * (1. + oii_coefficients(1:24,6) * (1. - tered) + oii_coefficients(1:24,7) * (1. - tered) ** 2))
      oii_coefficients(25:36,9)=1.e-14 * oii_coefficients(25:36,8) * tered**(oii_coefficients(25:36,5) + oii_coefficients(25:36,6)*log(tered) + oii_coefficients(25:36,7)* log(tered) ** 2)

! 4f-3d transitions

      if (tered .gt. 0.5) then
        aeff=oii_A_4f(9)
      else
        aeff=oii_A_4f_low(9)
      endif

      where (oiiRLs%Term1(4:5) .eq. "3d" .and. oiiRLs%Term2(3:4) .eq. "4f")
        oiiRLs%Em = aeff * 1.98648E-08 /oiiRLs%Wave * &
        & oiiRLs%g2 * oiiRLs%Br_B
        oiiRLs%Int = 100. * oiiRLs%Em / Em_hb * abund
      endwhere

! 3d-3p ^4F transitions (Case A=B=C for a,b,c,d; Br diff. slightly, adopt Case B)

      if (tered .gt. 0.5) then
        aeff=oii_A_3d4F(9)
      else
        aeff=oii_A_3d4F_low(9)
      endif

      where (oiiRLs%Mult .eq. " V10       " .or. oiiRLs%Mult .eq. " V18       ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_B
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3d-3p ^4D, ^4P transitions. case B assumed

      if (tered .gt. 0.5) then
        aeff=oii_B_3d4D(9)
      else
        aeff=oii_B_3d4D_low(9)
      endif

      where (oiiRLs%Term1(4:5) .eq. "3p" .and. (oiiRLs%Term2 .eq. "  3d  4D " .or. oiiRLs%Term2 .eq. "  3d  4P "))
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_B
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3d-3p ^2F transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3d2F(9)
      else
        aeff=oii_A_3d2F_low(9)
      endif

      where (oiiRLs%Term1(4:5) .eq. "3p" .and. oiiRLs%Term2 .eq. "  3d  2F   ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3d-3p ^2D transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3d2D(9)
      else
        aeff=oii_A_3d2D_low(9)
      endif

      where (oiiRLs%Term1(4:5) .eq. "3d" .and. oiiRLs%Term2 .eq. "  3d  2D   ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3d-3p ^2P transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3d2P(9)
      else
        aeff=oii_A_3d2P_low(9)
      endif

      where (oiiRLs%Term1(4:5) .eq. "3p" .and. oiiRLs%Term2 .eq. "  3d  2P   ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^4D - ^4P transitions. case B

      if (tered .gt. 0.5) then
        aeff=oii_B_3p4D(9)
      else
        aeff=oii_A_3p4D_low(9)
      endif

      where (oiiRLs%Mult .eq. " V1        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_B
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^4P - ^4P transitions. case B

      if (tered .gt. 0.5) then
        aeff=oii_B_3p4P(9)
      else
        aeff=oii_A_3p4P_low(9)
      endif
!
      where (oiiRLs%Mult .eq. " V2        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_B
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^4S - ^4P transitions. case B

      if (tered .gt. 0.5) then
        aeff=oii_B_3p4S(9)
      else
        aeff=oii_A_3p4S_low(9)
      endif
!
      where (oiiRLs%Mult .eq. " V3        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_B
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^2D - ^2P transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3p2D(9)
      else
        aeff=oii_A_3p2D_low(9)
      endif
!
      where (oiiRLs%Mult .eq. " V5        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^2P - ^2P transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3p2P(9)
      else
        aeff=oii_A_3p2P_low(9)
      endif
!
      where (oiiRLs%Mult .eq. " V6        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere
!
! 3p-3s ^2S - ^2P transitions. case A

      if (tered .gt. 0.5) then
        aeff=oii_A_3p2S(9)
      else
        aeff=oii_A_3p2S_low(9)
      endif
!
      where (oiiRLs%Mult .eq. " V4        ")
        oiiRLs%Em = aeff*1.98648E-08 / oiiRLs%Wave*&
      & oiiRLs%g2*oiiRLs%Br_A
        oiiRLs%Int = 100.*oiiRLs%Em / Em_hb*abund
      endwhere

end subroutine oii_rec_lines

subroutine read_oii_s2017

!read in tables of recombination coefficients from Storey et al. 2017.
!full file online at http://cdsarc.u-strasbg.fr/ftp/VI/150/DataFiles/OIIlines_ABC
!data for 276 lines extracted and put in Roii_new.txt
!data comes in blocks for log(Te)=2.0..0.1..4.4, and log(ne)=2.0..0.2..5.0

implicit none
real(kind=dp), dimension(17) :: temp
integer :: i,j,k,io
character(len=220) :: dump

!debugging
#ifdef CO
        print *,"subroutine: read_oii_s2017"
#endif

open(100,file=trim(PREFIX)//'/share/neat/Roii_storey2017.dat', iostat=IO, status='old')

do i=1,276
  read(100,"(A220)") dump ! first line of each block contains transition data

  read(dump(32:39),"(F7.2)") oiiRLs_s2017(i)%wave
  read(dump(129:134),"(A5)") oiiRLs_s2017(i)%mult
  read(dump(160:169),"(A9)") oiiRLs_s2017(i)%Term1
  read(dump(191:200),"(A9)") oiiRLs_s2017(i)%Term2

  read(100,*) dump ! second line is just restatement of density points

  do j=1,25 ! loop over temperatures
    read (100,*) temp ! read in row, 1 temperature value followed by 15 recombination coefficients
    do k=1,16 ! loop over density
      oiiRLs_s2017(i)%ems(j,k)=temp(k+1)
    enddo
  enddo

enddo

end subroutine read_oii_s2017

subroutine oii_rec_lines_s2017(Te,Ne,abund,oiiRLs_s2017,aeff_hb,em_hb)
! supersedes oii_rec_lines
! first interpolates logarithmically in Te and Ne to get all the recombination coefficients
! then calculates the predicted intensities

implicit none
real(kind=dp) :: aeff, aeff_hb, em_hb, Te, Ne, abund, tered
real(kind=dp) :: logte,logne
real(kind=dp) :: aa1, aa2, bb1, bb2, interp_te, interp_ne, em_t1, em_t2
integer :: te1,ne1,te2,ne2
integer :: i

type(oiiRL_s2017), dimension(276) :: oiiRLs_s2017

!debugging
#ifdef CO
        print *,"subroutine: oii_rec_lines_s2017"
#endif

! oii coefficients cover 16 densities (log(ne)=2.0..5.0..0.2) and 25 temperatures (log(te)=2.0..4.4..0.1)
! find the array indices to interpolate between
! 2.0 = 1    (log(te)-1.9)/0.1
! 2.1 = 2
! 2.2 = 3

te1=floor((log10(te)-1.9)/0.1)
ne1=floor((log10(ne)-1.8)/0.2)

if (te1 .lt. 1) then
  te1=1
  te2=1
  interp_te=0.
elseif (te1 .ge. 25) then
  te1=25
  te2=25
  interp_te=0.
else
  te2=te1+1
  interp_te=(log10(te)-1.9)/0.1 - real(te1)
endif

if (ne1 .lt. 1) then
  ne1=1
  ne2=1
elseif (ne1 .ge. 16) then
  ne1=16
  ne2=16
else
  ne2=ne1+1
  interp_ne=(log10(ne)-1.8)/0.2 - real(ne1)
endif

do i=1,size(oiiRLs_s2017)
  aa1=oiiRLs_s2017(i)%ems(te1,ne1)
  aa2=oiiRLs_s2017(i)%ems(te2,ne1)
  bb1=oiiRLs_s2017(i)%ems(te1,ne2)
  bb2=oiiRLs_s2017(i)%ems(te2,ne2)


  em_t1 = aa1+(aa2-aa1)*interp_te
  em_t2 = bb1+(bb2-bb1)*interp_te

  ! tabulated values are emission coefficients in erg.cm3/s
  ! need recombination coefficients in cm3/s
  ! = em / hc/lambda (in cgs) = 1.98644582e-25m3kg/lambda(m) = 1.986e-16cm3g/m = 1.986e-18cm3g/cm

  oiiRLs_s2017(i)%em_interpolated = em_t1+(em_t2-em_t1)*interp_ne
  oiiRLs_s2017(i)%Int = 100.*abund*oiiRLs_s2017(i)%em_interpolated / Em_hb
enddo

end subroutine oii_rec_lines_s2017

subroutine nii_rec_lines(te, ne, abund, niiRLs,aeff_hb,em_hb)

      implicit none
      real(kind=dp) :: aeff_hb, em_hb, Te, Ne, abund, tered

      type(niiRL), dimension(:) :: niiRLs

!debugging
#ifdef CO
        print *,"subroutine: nii_rec_lines"
#endif

      tered = te/10000.

      nii_coefficients(1:25,6) = 10. ** (nii_coefficients(1:25,1) + nii_coefficients(1:25,2) * log10(tered) + nii_coefficients(1:25,3) * log10(tered) ** 2)
      nii_coefficients(26:27,6) = (1.e-13 * 2.d0 * nii_coefficients(26:27,1) * (tered/4.d0)**nii_coefficients(26:27,2)) / (1.d0 + nii_coefficients(26:27,3)*(tered/4.d0)**nii_coefficients(26:27,4))*nii_coefficients(26:27,5)

!     2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p E1 3P* - 3D  M03  transitions
!     case B
!
      where (niiRLs%Mult .eq. "V3         ")
        niiRLs%Em = nii_V3(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!      2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3S     M04 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V4         ")
        niiRLs%Em = nii_V4(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 3P* - 3P     M05 transitions
!      case B
!

      where (niiRLs%Mult .eq. "V5         ")
        niiRLs%Em = nii_V5(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1P     M08 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V8         ")
        niiRLs%Em = nii_V8(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1D     M12 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V12        ")
        niiRLs%Em = nii_V12(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3s - 2s2.2p.(2P*).3p 1P* - 1S     M13 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V13        ")
        niiRLs%Em = nii_V13(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1D*     M15 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V15        ")
        niiRLs%Em = nii_V15(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1P - 1P*     M17 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V17        ")
        niiRLs%Em = nii_V17(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3F*     M19 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V19        ")
        niiRLs%Em = nii_V19(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3D*     M20 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V20        ")
        niiRLs%Em = nii_V20(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3D - 3P*     M21 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V21        ")
        niiRLs%Em = nii_V21(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3D - 3P*     M22 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V22        ")
        niiRLs%Em = nii_V22(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3S - 3P*     M24 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V24        ")
        niiRLs%Em = nii_V24(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3S - 3P*     M26 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V26        ")
        niiRLs%Em = nii_V26(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3D*     M28 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V28        ")
        niiRLs%Em = nii_V28(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 3P - 3P*     M29 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V29        ")
        niiRLs%Em = nii_V29(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).4s 3P - 3P*     M30 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V30        ")
        niiRLs%Em = nii_V30(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3p - 2s2.2p.(2P*).3d 1D - 1F*     M31 transitions
!      case A
!
      where (niiRLs%Mult .eq. "V31        ")
        niiRLs%Em = nii_V31(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*).4p 3F* - 3D     M36 transitions
!      case B
!
      where (niiRLs%Mult .eq. "V36        ")
        niiRLs%Em = nii_V36(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 3F* - 3G M39 transitions
!      case B
!
      where (niiRLs%Mult(1:3) .eq. "V39")
        niiRLs%Em = nii_V39(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       2s2.2p.(2P*).3d - 2s2.2p.(2P*<3/2>).4f 1F* - 1G M58 transitions
!      case A
!
      where (niiRLs%Mult(1:3) .eq. "V58")
        niiRLs%Em = nii_V58(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 3D* - 4f 3F 4242 M48 transitions
!      case B
!
      where (niiRLs%Mult(1:3) .eq. "V48")
        niiRLs%Em = nii_V48(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 3P* - 4f 3D 4435 M55 transitions
!      case B
!
      where (niiRLs%Mult(1:3) .eq. "V55")
        niiRLs%Em = nii_V55(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 1D* - 4f 1F 4176 M43 (RMT M42) transitions
!      case A
!
      where (niiRLs%Mult(1:3) .eq. "V43")
        niiRLs%Em = nii_V43(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 1P* - 4f 1D 4677 M61 (RMT M62) transitions
!      case A
!
      where (niiRLs%Mult(1:3) .eq. "V61")
        niiRLs%Em = nii_V61(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 3F* - 4f 1G 4026 M39b transitions
!      case A (PPB):
!
      where (niiRLs%Mult(1:4) .eq. "V39b")
        niiRLs%Em = nii_V39b(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere
!
!       3d 1F* - 4f 3G 4552 M58a transitions
!      case A (PPB):
!
      where (niiRLs%Mult(1:4) .eq. "V58a")
        niiRLs%Em = nii_V58a(6) * 1.98648E-08 / niiRLs%Wave * niiRLs%Br_LS
        niiRLs%Int = 100 * niiRLs%Em / Em_Hb * abund
      endwhere

end subroutine nii_rec_lines

subroutine cii_rec_lines(te, ne, abund, ciiRLs,aeff_hb,em_hb)

      implicit none
      real(kind=dp) :: aeff_Hb, em_hb, Te, Ne, abund, tered

      type(ciiRL), dimension(:) :: ciiRLs

!debugging
#ifdef CO
        print *,"subroutine: cii_rec_lines"
#endif

      tered = te/10000

      ciiRLs%aeff = 1e-14 * (ciiRLs%a*(tered**ciiRLs%f)) * (1 &
      &+ (ciiRLs%b*(1-tered)) &
      &+ (ciiRLs%c * ((1-tered)**2) ) &
      &+ (ciiRLs%d * ((1-tered)**3) ) &
      &)
      ciiRLs%Int = 100 * (ciiRLs%aeff/aeff_hb) * (4861.33/ciiRLs%Wave) * abund

end subroutine cii_rec_lines

subroutine neii_rec_lines(te, ne, abund, neiiRLs,aeff_hb,em_hb)

      implicit none
      real(kind=dp) :: aeff_Hb, em_hb, Te, Ne, abund, tered

      type(neiiRL), dimension(:) :: neiiRLs

!debugging
#ifdef CO
        print *,"subroutine: neii_rec_lines"
#endif

      tered = te/10000.

      neiiRLs%aeff = neiiRLs%Br * 1e-14 * &
      &(neiiRLs%a*(tered**neiiRLs%f)) * (1 &
      &+ (neiiRLs%b*(1-tered)) &
      &+ (neiiRLs%c * ((1-tered)**2) ) &
      &+ (neiiRLs%d * ((1-tered)**3) ) &
      &)
      neiiRLs%Int = 100 * (neiiRLs%aeff/aeff_hb) * (4861.33/neiiRLs%Wave) * abund

end subroutine neii_rec_lines

subroutine xiii_rec_lines(te, ne, abund, xiiiRLs,aeff_hb,em_hb)

      implicit none
      real(kind=dp) :: aeff_Hb, em_hb, Te, Ne, abund, tered

      type(xiiiRL), dimension(:) :: xiiiRLs

      tered = te/90000. !ionic charge=3 so divide by 9

      xiiiRLs%aeff = xiiiRLs%Br * 1e-13 * 3 * &
      & (xiiiRLs%a*(tered**xiiiRLs%b)) / &
      & (1 + (xiiiRLs%c * (tered**xiiiRLs%d)))
      xiiiRLs%Int = 100 * (xiiiRLs%aeff/aeff_hb) * (4861.33/xiiiRLs%Wave) * abund

end subroutine xiii_rec_lines

subroutine get_aeff_hb(te, ne, aeff_hb, em_hb)
      implicit none
      real(kind=dp) :: Te, Ne, AE2, AE3, AE4, AE5, AE6, AE7, AE8, aeff_hb, Em_Hb, logem

!debugging
#ifdef CO
        print *,"subroutine: get_aeff_hb"
#endif

      AE2 = -9.06524E+00 -2.69954E+00 * log10(te) + 8.80123E-01 * &
      &log10(te) ** 2 -1.57946E-01 * log10(te) ** 3 + &
      &9.25920E-03 * log10(te) ** 4
      AE3 = -8.13757E+00 -3.57392E+00 * log10(te) + 1.19331E+00 * &
      &log10(te) ** 2 -2.08362E-01 * log10(te) ** 3 + &
      &1.23303E-02 * log10(te) ** 4
      AE4 = -6.87230E+00 -4.72312E+00 * log10(te) + 1.58890E+00 * &
      &log10(te) ** 2 -2.69447E-01 * log10(te) ** 3 + &
      &1.58955E-02 * log10(te) ** 4
      AE5 = -5.15059E+00 -6.24549E+00 * log10(te) + 2.09801E+00 * &
      &log10(te) ** 2 -3.45649E-01 * log10(te) ** 3 + &
      &2.01962E-02 * log10(te) ** 4
      AE6 = -2.35923E+00 -8.75565E+00 * log10(te) + 2.95600E+00 * &
      &log10(te) ** 2 -4.77584E-01 * log10(te) ** 3 + &
      &2.78852E-02 * log10(te) ** 4
      AE7 =  1.55373E+00 -1.21894E+01 * log10(te) + 4.10096E+00 * &
      &log10(te) ** 2 -6.49318E-01 * log10(te) ** 3 + &
      &3.76487E-02 * log10(te) ** 4
      AE8 =  6.59883E+00 -1.64030E+01 * log10(te) + 5.43844E+00 * &
      &log10(te) ** 2 -8.40253E-01 * log10(te) ** 3 + &
      &4.79786E-02 * log10(te) ** 4

      if (log10(ne) .lt. 2) then
            aeff_hb = ae2
      elseif (log10(ne) .ge. 2 .and. log10(ne) .lt. 3) then
            aeff_hb = AE2 + (AE3 - AE2) * (log10(ne) - 2)
      elseif (log10(ne) .ge. 3 .and. log10(ne) .lt. 4) then
            aeff_hb = AE3 + (AE4 - AE3) * (log10(ne) - 3)
      elseif (log10(ne) .ge. 4 .and. log10(ne) .lt. 5) then
            aeff_hb = AE4 + (AE5 - AE4) * (log10(ne) - 4)
      elseif (log10(ne) .ge. 5 .and. log10(ne) .lt. 6) then
            aeff_hb = AE5 + (AE6 - AE5) * (log10(ne) - 5)
      elseif (log10(ne) .ge. 6 .and. log10(ne) .lt. 7) then
            aeff_hb = AE6 + (AE7 - AE6) * (log10(ne) - 6)
      elseif (log10(ne) .ge. 7 .and. log10(ne) .lt. 8) then
            aeff_hb = AE7 + (AE8 - AE7) * (log10(ne) - 7)
      else
            aeff_hb = AE8
      endif

      LogEm = aeff_hb - 11.38871 ! = log10(hc/lambda in cgs)
      aeff_hb = 10**aeff_hb
      em_hb = 10**logem

end subroutine get_aeff_hb

end module mod_recombination_lines
