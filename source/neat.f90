! NEAT, the nebular abundance analysis tool
! (C) 2006- Roger Wesson, Dave Stock, Peter Scicluna
! NEAT incorporates aspects of several codes developed over decades at UCL:
! equib, by I. Howarth and S. Adams, updated significantly by B. Ercolano
! MIDAS scripts written by X-W. Liu for calculating recombination line abundances
! It also uses the quicksort algorithm as implemented in F90 by Alberto Ramos

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

program neat

        use mod_types
        use mod_globals
        use mod_resultarrays
        use mod_extinction
        use mod_quicksort
        use mod_abundIO
        use mod_atomicdata
        use mod_atomic_read
        use mod_recombination_lines
        use mod_linefinder
        use mod_hydrogen
        use mod_helium
        use mod_oii_diagnostics
        use mod_weights
        use mod_functions
        use mod_commandline
        use mod_output

        implicit none

        character :: switch_ext !switch for extinction laws
        character :: switch_he  !switch for helium atomic data
        character :: switch_icf !switch for which ICF scheme to use
        integer :: I, j, runs, Narg !runs = number of runs for randomiser

!file reading variables

        logical :: identifylines,identifyconfirm
        type(line),dimension(:), allocatable :: linelist
        type(line),dimension(:), allocatable :: linelist_original
        type(line),dimension(:,:), allocatable :: all_linelists
        character(len=1) :: blank
        integer :: listlength, ncols, strpos
        real(kind=dp) :: normalise
        logical :: fitsinput

!results and result processing

        type(resultarray), dimension(:), allocatable :: all_results
        type(resultarray), dimension(1) :: iteration_result
        real(kind=dp), dimension(:), allocatable :: quantity_result

!atomic data

        character(len=10), dimension(24) :: ionlist !list of ion names
        integer :: iion !# of ions in Ilines
        integer :: maxlevs,maxtemps
        type(atomic_data),allocatable :: atomicdata(:)
        real(kind=dp), dimension(:,:,:), allocatable :: heidata

!line indexing arrays - these contain the location within the main array of the lines of interest, plus levels and zones for CELs

        integer, dimension(3:40) :: H_Balmer
        integer, dimension(4:39) :: H_Paschen
        integer, dimension(44) :: HeI_lines
        integer, dimension(20,2:6) :: HeII_lines
        type(cel), dimension(88) :: ILs !todo: work out why this becomes undefined on entry if its shape is assumed

!extinction

        logical :: calculate_extinction
        real(kind=dp) :: meanextinction, R

!RP effect switch

        logical :: norp

!subtract recombination contribution to auroral lines switch - 1=no subtraction, 2=subtraction

        integer :: subtract_recombination

!diagnostics

        type(weightingarray) :: weights

!binning and uncertainties

        integer :: verbosity,nbins,nperbin

!OpenMP

        integer :: omp_get_num_threads

!diagnostic array

        type(diagnostic_array) :: diagnostics

! default values

        runs=1
        switch_ext="S" !Howarth 1983 Galactic law
        switch_he="P"  !Porter 2012 data
        switch_icf="D" !DI14 ICF
        filename=""
        meanextinction=0.D0
        diagnostics%lowtemp=0.D0
        diagnostics%lowdens=0.D0
        diagnostics%medtemp=0.D0
        diagnostics%meddens=0.D0
        diagnostics%hightemp=0.D0
        diagnostics%highdens=0.D0
        verbosity=3
        R=3.1
        identifylines=.false.
        identifyconfirm=.false.
        nbins=25
        normalise = 0.d0
        norp=.true.
        calculate_extinction=.true.
        subtract_recombination=1
        configfile=""
        defaultconfigfile=trim(PREFIX)//"/share/neat/default.cfg"
        fitsinput=.false.

        iion=24 ! number of ions for which we can calculate abundances
                ! todo: calculate this automatically
                ! and replace it in calls to subroutines with getting the array size within the subroutine

        print *,"NEAT, the Nebular Empirical Analysis Tool"
        print *,"version ",VERSION
        print *
        print *,gettime(),"starting code"

        call readcommandline(runs,switch_ext,switch_he,switch_icf,meanextinction,diagnostics,verbosity,R,identifylines,identifyconfirm,nbins,normalise,norp,calculate_extinction,subtract_recombination,configfile,nperbin)

! read in the line list, allocate original linelist array

        strpos=index(filename,".",.true.)+1
        if (trim(filename(strpos:len(filename))).eq."fit" .or. trim(filename(strpos:len(filename))).eq."fits".or.trim(filename(strpos:len(filename))).eq."FIT".or.trim(filename(strpos:len(filename))).eq."FITS") then
          fitsinput=.true.
          call read_fits_linelist(linelist,listlength,ncols,runs)
        else
          call read_text_linelist(linelist,listlength,ncols,runs)
        endif

        allocate (linelist_original(listlength))

! run line identifier if required

        if (identifylines) then
                call linefinder(linelist, listlength, identifyconfirm)
        endif

! normalise to Hbeta

        if (minval(abs(linelist(:)%wavelength - 4861.33)) .lt. 0.001) then
          normalise = 100.d0/linelist(minloc(abs(linelist(:)%wavelength - 4861.33),1))%intensity
        else
          print *,gettime(),"error: no H beta detected. No further analysis possible"
          call exit(1)
        endif

        linelist%intensity = linelist%intensity * normalise
        linelist%int_err = linelist%int_err * normalise
        linelist%blend_intensity = linelist%blend_intensity * normalise
        linelist%blend_int_err = linelist%blend_int_err * normalise

! check for and remove negative line fluxes and uncertainties

        if (minval(linelist%intensity) .lt. 0.) then
          where (linelist%intensity .lt. 0)
            linelist%intensity = 0.D0
          endwhere
        print *,gettime(),"warning: negative line fluxes set to zero"
        endif

        if (minval(linelist%int_err) .lt. 0.) then
          where (linelist%int_err .lt. 0)
            linelist%int_err = abs(linelist%int_err)
          endwhere
        print *,gettime(),"warning: negative flux uncertainties converted to positive"
        endif

        print *

! select reddening law, get flambda from functions if one of the five, otherwise read from file

        if (switch_ext == "S") then
                print *,gettime(),"using Howarth (1983) galactic extinction law"
        elseif (switch_ext == "H") then
                print *,gettime(),"using Howarth (1983) LMC extinction law"
        elseif (switch_ext == "C") then
                print *,gettime(),"using CCM (1989) galactic extinction law"
        elseif (switch_ext == "P") then
                print *,gettime(),"using Prevot et al. (1984) SMC extinction law"
        elseif (switch_ext == "F") then
                print *,gettime(),"using Fitzpatrick (1990) galactic extinction law"
        endif

        call get_flambda(linelist, switch_ext, R)

        if (switch_he == "S") then
                print *,gettime(),"using Smits (1996) He I emissivities"
        else
                print *,gettime(),"using Porter et al. (2012) He I emissivities"
        endif

        if (switch_ICF == "K") then
                print *,gettime(),"using Kingsburgh & Barlow (1994) ICF"
        elseif (switch_ICF == "D") then
                print *,gettime(),"using Delgado-Inglada et al. (2014) ICF"
        elseif (switch_ICF == "E") then
                print *,gettime(),"using Delgado-Inglada et al. (2014) ICF with classical icf(N)"
        else
                print *,gettime(),"using Peimbert et al. (1992) ICF"
        endif

        print *

! read the CEL data

        call read_celdata(ILs,ionlist)

! do line identifying here

        !first find H & He lines, and CELs
        call get_H(H_Balmer, H_Paschen, linelist)
        call get_Hei(Hei_lines, linelist)
        call get_Heii(Heii_lines, linelist)
        call get_cels(ILs,linelist)

! read in the weights to be used in the analysis.

        print *,gettime(),"reading in abundance analysis weights from ",trim(defaultconfigfile)
        call setweights(defaultconfigfile,weights,linelist,ILs,H_Balmer,H_Paschen,HeII_lines)
        if (configfile .ne. "") then
          print *,gettime(),"reading in abundance analysis weights from ",trim(configfile)
          print *,gettime(),"(defaults still apply for any weights not specified in this file)"
          call setweights(configfile,weights,linelist,ILs,H_Balmer,H_Paschen,HeII_lines)
        endif

!CELs read, can now read in atomic data

        print *,gettime(),"reading atomic data from ",trim(PREFIX),"/share/neat"
        allocate(atomicdata(iion))
        do I=1,iion
            atomicdata(I)%ion = ionlist(I)
            call read_atomic_data(atomicdata(I))
        enddo

!read ORL data
        call read_orl_data

!read hydrogen emissivities

        print *,gettime(),"reading H emissivities"
        call read_hydrogen

!read helium emissivities

        print *,gettime(),"reading He I emissivities"
        if (switch_he .eq. "P") then
          allocate(heidata(21,14,44))
          heidata = 0.D0
          call read_porter(heidata)
        elseif (switch_he .eq. "S") then
          allocate(heidata(3,6,44))
          heidata = 0.D0
          call read_smits(heidata)
        endif

        print *,gettime(),"reading He II emissivities"
        call read_heii

!read oii data for diagnostics

        print *,gettime(),"reading OII diagnostics data"
        call read_oii_diagnosticdata
        call read_oii_s2017

!find maximum #levels and temperatures - pass to equib to reduce footprint

        maxlevs = maxval(atomicdata%nlevs)
        maxtemps = maxval(atomicdata%ntemps)

!now check number of iterations.  If 1, line list is fine as is.  If more than one, randomize the fluxes

        allocate(all_results(runs))
        allocate(all_linelists(size(linelist),runs))

        if (runs == 1)then !calculates abundances without uncertainties
                print *
                print *,gettime(),"doing abundance calculations"

                call abundances(linelist, listlength, iteration_result, meanextinction, calculate_extinction, ILs, diagnostics,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf, H_Balmer, H_Paschen, HeI_lines, HeII_lines, weights,subtract_recombination)

                all_results(1)=iteration_result(1) ! copy the iteration result to all_results to simplify the writing out of results later
                all_linelists(:,1)=linelist

                print *,gettime(),"finished abundance calculations"
        else if (runs .gt. 1)then

!save unrandomised line list

                linelist_original = linelist

                call init_random_seed()!sets seed for randomiser

                !main loop
!$OMP PARALLEL default(firstprivate) shared(all_linelists,listlength,norp,calculate_extinction,ILs,diagnostics,iion,atomicdata,maxlevs,maxtemps,switch_he,switch_icf,all_results,subtract_recombination)
!$OMP MASTER

                print *
                if (omp_get_num_threads().gt.1) then
                  print "(X,A10,X,A,I2,A)",gettime(),"starting Monte Carlo calculations, using ",omp_get_num_threads()," processors"
                else
                  print "(X,A10,X,A)",gettime(),"starting Monte Carlo calculations"
                endif
                print *,gettime(),"completed ",0,"%"
!$OMP END MASTER
!$OMP DO schedule(dynamic)
                do I=1,runs

                         if ( (10.0*dble(i)/dble(runs)) == int(10*i/runs) ) print *,gettime(),"completed ",100*i/runs,"%"
!                        print*, "iteration ", i, "of", runs

                        call randomizer(linelist, listlength, norp)
                        call abundances(linelist, listlength, iteration_result, meanextinction, calculate_extinction, ILs, diagnostics,iion,atomicdata,maxlevs,maxtemps, heidata, switch_he, switch_icf, H_Balmer, H_Paschen, HeI_lines, HeII_lines, weights, subtract_recombination)

                        !store all line and derived quantity in arrays
                        all_linelists(:,i)=linelist
                        all_results(i)=iteration_result(1)
                        !restore the unrandomized line list ready for the next iteration
                        linelist = linelist_original
                enddo
!$OMP END DO
!$OMP END PARALLEL
        else
                print*, gettime(),"error: I didn't want to be a barber anyway. I wanted to be... a lumberjack!   Also, a positive number of runs helps.."
                call exit(1)
        endif

!now write out the line list and results files

        if (fitsinput) then
          call write_fits(runs,listlength,ncols,all_linelists,all_results,verbosity,nbins,subtract_recombination)
        else
          call write_output(runs,listlength,ncols,all_linelists,all_results,verbosity,nbins,subtract_recombination)
        endif

        print *
        print *,gettime(),"all done"
        print *

end program
