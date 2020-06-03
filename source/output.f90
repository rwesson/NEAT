!output.f90, for writing out the linelist and results
!(C) Roger Wesson
module mod_output
use mod_functions
use mod_globals

implicit none

contains

subroutine write_output(runs,listlength,ncols,all_linelists,all_results,verbosity,nbins,subtract_recombination)

        type(line), dimension(:,:) :: all_linelists
        type(resultarray), dimension(:) :: all_results
        type(resultarray), dimension(1) :: iteration_result
        real(kind=dp), dimension(:), allocatable :: quantity_result
        real(kind=dp), dimension(:,:), allocatable :: resultprocessingarray
        character(len=40), dimension(:,:), allocatable :: resultprocessingtext
        integer :: runs,listlength,ncols,verbosity
        real(kind=dp), dimension(3) :: uncertainty_array=0d0
        type(arraycount), dimension (:), allocatable :: binned_quantity_result
        logical :: unusual
        integer :: nbins,subtract_recombination
        integer :: i,j

#ifdef CO
        print *,"subroutine: write_output"
#endif

        print *
        print *,gettime(),"Writing line list"

        allocate(quantity_result(runs))
        quantity_result=0d0
        open (650,file=trim(filename)//"_linelist", status='replace', access='sequential', action='write')
        open (651,file=trim(filename)//"_linelist.tex", status='replace', access='sequential', action='write')

        if (ncols .ge. 4) then
            write (650,*) " Lambda    (Rest)  Ion         F(line)  I(line)      Abundance"
            write (651,*) "\begin{longtable}{llrrlllllll}"
            write (651,*) "\hline"
            write (651,*) "$\lambda_{obs}$ & $\lambda_{rest}$ & $F \left( \lambda \right) $ & $I \left( \lambda \right) $ & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
        else
            write (650,*) " Lambda  Ion         F(line)  I(line)      Abundance"
            write (651,*) "\begin{longtable}{lrrlllllll}"
            write (651,*) "\hline"
            write (651,*)                   "$\lambda_{rest}$ & $F \left( \lambda \right) $ & $I \left( \lambda \right) $ & Ion & Multiplet & Lower term & Upper term & g$_1$ & g$_2$ \\"
        endif
        write (651,*) "\hline"

        if (runs .gt. 1) then
                do j=1, listlength

!observed wavelength if known

                if (ncols .ge. 4) then
                  if (all_linelists(j,1)%wavelength_observed .gt. 0.) then
                    write (650,"(X,F8.2,X)", advance='no') all_linelists(j,1)%wavelength_observed
                    write (651,"(X,F8.2,' & ')", advance='no') all_linelists(j,1)%wavelength_observed
                  else
                    write (650,"(X,A,X)", advance='no') "       *"
                    write (651,"(X,A,' & ')", advance='no') "       *"
                  endif
                endif

!rest wavelength, ion name for plain text file

                write (650,"(X,F8.2,X,A11)", advance='no') all_linelists(j,1)%wavelength,all_linelists(j,1)%name
                write (651,"(X,F8.2,' & ',A15,' & ')", advance='no') all_linelists(j,1)%wavelength

!line flux

               if (all_linelists(j,1)%intensity .eq. 0.d0 .and. all_linelists(j,1)%blend_intensity .eq. 0.d0) then
                 write (650,"(A)", advance='no') " * "
                 write (651,"(A)", advance='no') " *     &             &"
               else
                 write (650,"(F8.3,A,F8.3,3X)", advance='no') all_linelists(j,1)%intensity," +-",all_linelists(j,1)%int_err
                 write (651,"(F8.3,'& $\pm$',F8.3, '&')", advance='no') all_linelists(j,1)%intensity,all_linelists(j,1)%int_err
               endif

!dereddened flux
                quantity_result = all_linelists(j,:)%int_dered

                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)

                if (all_linelists(j,1)%intensity .ne. 0.d0 .or. all_linelists(j,1)%blend_intensity .ne. 0.d0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(F8.3,SP,F8.3,SP,F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                    write (651,"(F8.3,'& $^{',SP,F8.3,'}_{',SP,F8.3,'}$')", advance='no') uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
                  else
                    write (650,"(F8.3,A,F8.3,4X)", advance='no') uncertainty_array(2)," +-",uncertainty_array(1)
                    write (651,"(F8.3,'& $\pm$',F8.3)", advance='no') uncertainty_array(2),uncertainty_array(1)
                  endif
                else
                  write (650,"(A)", advance='no') " * "
                  write (651,"(A)", advance='no') " *     &        "
                endif

! transition data

!                write (650,*)
                write (651,*) all_linelists(j,1)%linedata, "\\"

!abundance - write out if there is an abundance for the line, don't write
!anything except a line break if there is no abundance for the line.
!todo: add an option to choose whether or not to put abundances in the line list table

                quantity_result = all_linelists(j,:)%abundance
                call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)
                if (uncertainty_array(2) .ne. 0.D0) then
                  if (uncertainty_array(1) .ne. uncertainty_array(3)) then
                    write (650,"(ES10.2,SP,ES10.2,SP,ES10.2)") uncertainty_array(2),uncertainty_array(1),-uncertainty_array(3)
!                    write (651,"(' & ${',A,'}$ & $^{+',A,'}_{',A,'}$ \\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1))),trim(latex_number(-uncertainty_array(3)))
                  else
                    write (650,"(ES10.2,A,ES10.2)") uncertainty_array(2)," +-",uncertainty_array(1)
!                    write (651,"(' & $',A,'$ & $\pm',A,'$\\')") trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
                  endif
                else
                  write (650,*)
!                  write (651,*) "\\"
                endif

                enddo
        else ! runs == 1, no uncertainties to write out

                do i=1,listlength
                  if (ncols .ge. 4) then
                    if (all_linelists(i,1)%wavelength_observed .gt. 0.) then
                      write (650,"(X,F8.2,X)", advance='no') all_linelists(i,1)%wavelength_observed
                      write (651,"(X,F8.2,' & ')", advance='no') all_linelists(i,1)%wavelength_observed
                    else
                      write (650,"(X,A,X)", advance='no') "       *"
                      write (651,"(X,A,' & ')", advance='no') "       *"
                    endif
                  endif
                  if (all_linelists(i,1)%intensity .ne. 0.d0) then
                    if (all_linelists(i,1)%abundance .gt. 0.0) then
                      write (650,"(X,F8.2,X,A11,F8.3,X,F8.3,X,ES14.3)") all_linelists(i,1)%wavelength,all_linelists(i,1)%name,all_linelists(i,1)%intensity,all_linelists(i,1)%int_dered, all_linelists(i,1)%abundance
                      write (651,"(X,F8.2,X,'&',X,F8.3,X,'&',X,F8.3,X,A,'\\')") all_linelists(i,1)%wavelength,all_linelists(i,1)%intensity,all_linelists(i,1)%int_dered,all_linelists(i,1)%linedata
                    else
                      write (650,"(X,F8.2,X,A11,F8.3,X,F8.3)") all_linelists(i,1)%wavelength,all_linelists(i,1)%name,all_linelists(i,1)%intensity,all_linelists(i,1)%int_dered
                      write (651,"(X,F8.2,X,'&',X,F8.3,X,'&',X,F8.3,X,A,'\\')") all_linelists(i,1)%wavelength,all_linelists(i,1)%intensity,all_linelists(i,1)%int_dered,all_linelists(i,1)%linedata
                    endif
                  else
                    write (650,"(X,F8.2,X,A11,A7,X,A7)") all_linelists(i,1)%wavelength,all_linelists(i,1)%name,"*      ","*      "
                    write (651,"(X,F8.2,X,'&',X,A7,X,'&',X,A7,X,A,'\\')") all_linelists(i,1)%wavelength,"*      ","*      ",all_linelists(i,1)%linedata
                  endif
                enddo

        endif

        write (651,*) "\hline"
        !write (651,*) "\caption{}"
        write (651,*) "\label{tab:",trim(filename)//"_linelist}"
        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),"Linelist written to files:"
        print *,"              ",trim(filename),"_linelist"
        print *,"              ",trim(filename),"_linelist.tex"

!now write out the summary files and all the binned data

        if (verbosity .eq. 1) then
          print *,gettime(),"Writing summary files, binned results and complete results"
        elseif (verbosity .eq. 2) then
         print *,gettime(),"Writing summary files and binned results"
        else
          print *,gettime(),"Writing summary files"
        endif

!get the arrays of output quantities

        call create_output_arrays(all_results,runs,resultprocessingarray,resultprocessingtext)

!open the files and write the headers

        open (650,file=trim(filename)//"_results", status='replace', access='sequential', action='write')
        open (651,file=trim(filename)//"_results.tex", status='replace', access='sequential', action='write')

        write (650,*) "NEAT (nebular empirical analysis tool)"
        write (650,*) "======================================"
        write (650,*) "Version ",VERSION
        write (650,*)
        write (650,*) "Analysis of file ",trim(filename)
        write (650,*) "Command line: ",trim(commandline)
        write (650,*)

        write (651,*) "\noindent{\Large {\sc neat} (nebular empirical analysis tool)}"
        write (651,*) "\noindent Version ",VERSION
        write (651,*) "\hrule"
        write (651,*) "\vspace{0.3cm}"
        write (651,*) "\noindent Analysis of file {\tt ",trim(filename),"}\newline"
        write (651,*) "\noindent Command line: {\tt ",trim(commandline),"}\newline"
        write (651,*) "\begin{longtable}[l]{ll}"

!next, loop through the results, processing and printing

        do j=1,164

! here we put some if statements to put things into conveniently separate bits

          if (j .eq. 1) then
            write (650,"(/A,/A/)") "Extinction","=========="
            write (651,*) "\multicolumn{2}{l}{Extinction}\\ \hline"
          elseif (j .eq. 5) then
            write (650,"(/A,/A)") "Diagnostics","==========="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Diagnostics}\\ \hline"
            write (650,"(/A,/A/)") "Low ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation densities}\\ \hline"
          elseif (j .eq. 10) then
            write (650,"(/A,/A/)") "Low ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Low ionisation temperatures}\\ \hline"
          elseif (j .eq. 21) then
            write (650,"(/A,/A/)") "Medium ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation densities}\\ \hline"
          elseif (j .eq. 36) then
            write (650,"(/A,/A/)") "Medium ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Medium ionisation temperatures}\\ \hline"
          elseif (j .eq. 51) then
            write (650,"(/A,/A/)") "High ionisation densities","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation densities}\\ \hline"
          elseif (j .eq. 54) then
            write (650,"(/A,/A/)") "High ionisation temperatures","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{High ionisation temperatures}\\ \hline"
          elseif (j .eq. 59) then
            write (650,"(/A,/A/)") "Recombination line diagnostics","-----------"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination line diagnostics}\\ \hline"
          elseif (j .eq. 71) then
            if (subtract_recombination .eq. 2) then
              write (650,"(/A,/A/)") "Recombination contribution to CELs (%) (RL value has been subtracted from line fluxes)","-----------"
              write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination contribution to CELs (\%)}\\ \hline"
            else
              write (650,"(/A,/A/)") "Recombination contribution to CELs (%) (for information only, not subtracted from line fluxes)","-----------"
              write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Recombination contribution to CELs (\%)}\\ \hline"
            endif
          elseif (j .eq. 77) then
            write (650,"(/A,/A/)") "CEL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{CEL abundances}\\ \hline"
          elseif (j .eq. 115) then
            write (650,"(/A,/A/)") "ORL abundances","=============="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{ORL abundances}\\ \hline"
          elseif (j .eq. 83 .or. j .eq. 89 .or. j .eq. 95 .or. j .eq. 101 .or. j .eq. 106 .or. j .eq. 110 .or. j .eq. 115 .or. j .eq. 118 .or. j .eq. 122 .or. j .eq. 133 .or. j .eq. 148) then
            if (j .eq. 133 .and. any(all_results%niiRLreliable .eqv. .false.)) then
              write (650,*) "(multiplet abundances do not agree - abundance may be unreliable)"
            endif
            if (j .eq. 148 .and. any(all_results%oiiRLreliable .eqv. .false.)) then
              write (650,*) "(multiplet abundances do not agree and/or V1 and V10 not detected - abundance may be unreliable)"
            endif
            write (650,"(/A,/A/)")
            write (651,*) "\\"
          elseif (j .eq. 151) then
            write (650,"(/A,/A/)") "Strong line abundances","======================"
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Strong line abundances}\\ \hline"
          elseif (j .eq. 157) then
            write (650,"(/A,/A/)") "Abundance discrepancy factors","============================="
            write (651,*) "\vspace{0.2cm}\\\multicolumn{2}{l}{Abundance discrepancy factors}\\ \hline"
          endif

! this writes the results to the plain text and latex summary files

          quantity_result=resultprocessingarray(j,:)
          call write_uncertainties(quantity_result,uncertainty_array,resultprocessingtext(j,1),resultprocessingtext(j,2),resultprocessingtext(j,3),filename, resultprocessingtext(j,4), verbosity,nbins)

        enddo

!write ends of files, close

        write (651,*) "\end{longtable}"

        close(650)
        close(651)

        print *,gettime(),"Results written to files:"
        print *,"              ",trim(filename),"_results"
        print *,"              ",trim(filename),"_results.tex"

end subroutine write_output

subroutine write_fits(runs,listlength,ncols,all_linelists,all_results,verbosity,nbins,subtract_recombination)
! fits output. requires input to have been ALFA-generated FITS, which will have four columns.

  implicit none
  type(line), dimension(:,:) :: all_linelists
  type(resultarray), dimension(:) :: all_results
  type(resultarray), dimension(1) :: iteration_result
  real(kind=dp), dimension(:), allocatable :: quantity_result
  integer :: runs,listlength,ncols,verbosity
  integer :: nbins,subtract_recombination
  real(kind=dp), dimension(:,:), allocatable :: resultprocessingarray
  character(len=40), dimension(:,:), allocatable :: resultprocessingtext
  integer :: i,j
  real(kind=dp), dimension(3) :: uncertainty_array=0d0
  type(arraycount), dimension (:), allocatable :: binned_quantity_result
  logical :: unusual
  character(len=30) :: cfitsioerror

!cfitsio variables
  integer :: status,unit,readwrite,blocksize,tfields,varidat
  character(len=16) :: extname
  character(len=16),dimension(6) :: ttype_lines,tform_lines,tunit_lines
  character(len=16),dimension(4) :: ttype_results,tform_results,tunit_results

#ifdef CO
        print *,"subroutine: write_fits"
#endif

  status=0
  readwrite=1

! open the input file, add date and history to primary HDU

  call ftgiou(unit,status)
  call ftopen(unit,trim(filename),readwrite,blocksize,status)
  call ftpdat(unit,status)
  call ftphis(unit,trim(commandline),status)

! go to extension LINES to update linelist

  call ftmnhd(unit,-1,"LINES",0,status)

! todo: check number of columns, overwrite if previous analysis exists

! add columns: dereddened flux, err, ion, linedata,abundance,err

  ttype_lines=(/"DereddenedFlux  ","DereddenedFluxE ","Ion             ","Linedata        ","Abundance       ","AbundanceUncert "/)
  tform_lines=(/"1E ","1E ","10A","85A","1E ","1E "/)
  tunit_lines=(/"                ","                ","                ","                ","                ","                "/)

  call fticls(unit,7,6,ttype_lines,tform_lines,status)
  call ftpcls(unit,9,1,1,listlength,all_linelists(:,1)%name,status)
  call ftpcls(unit,10,1,1,listlength,all_linelists(:,1)%linedata,status)

! loop to get uncertainties
! todo: maybe quicker to add property to type in the loop, then write to FITS in one go afterwards?
! allocate quantity result?
  do i=1,listlength
    quantity_result = all_linelists(i,:)%int_dered
    call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)
    call ftpcld(unit,7,i,1,1,uncertainty_array(2),status)

    quantity_result = all_linelists(i,:)%abundance
    call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)
    call ftpcld(unit,11,i,1,1,uncertainty_array(2),status)
  enddo

  print *,gettime(),"updated LINES extension"

! if RESULTS extension does not exist, create it, otherwise overwrite
! columns for quantity name, value, upper and lower uncertainties
! 162 quantities calculated

  tfields=4
  extname="RESULTS"
  ttype_results=(/"Quantity        ","Value           ","UpperUncertainty","LowerUncertainty"/)
  tform_results=(/"40A","1E ","1E ","1E "/)
  tunit_results=(/"                ","                ","                ","                "/)

  call ftmnhd(unit,-1,"RESULTS",0,status)
  if (status.eq.301) then
    print *,gettime(),"created RESULTS extension"
    status=0
    call ftmnhd(unit,-1,"QC",0,status)
    call ftibin(unit,164,tfields,ttype_results,tform_results,tunit_results,extname,varidat,status)
  else
    print *,gettime(),"RESULTS extension already present - overwriting"
  endif

! add header comments

  call ftpcom(unit,"Produced by neat version "//VERSION,status)
  call ftpcom(unit,"Command line: '"//trim(commandline)//"'",status)
  call ftpcom(unit,"input file: "//trim(filename),status)

! todo: record extinction law, helium data, ICF

! get the result arrays

  call create_output_arrays(all_results,runs,resultprocessingarray,resultprocessingtext)

! write out
! get_uncertainties returns array with (-,value,+)

  do j=1,164
    quantity_result=resultprocessingarray(j,:)
    call get_uncertainties(quantity_result, binned_quantity_result, uncertainty_array, unusual,nbins)
    call ftpcls(unit,1,j,1,1,resultprocessingtext(j,1),status)
    call ftpcld(unit,2,j,1,1,uncertainty_array(2),status)
    call ftpcld(unit,3,j,1,1,uncertainty_array(2)+uncertainty_array(1),status)
    call ftpcld(unit,4,j,1,1,uncertainty_array(2)-uncertainty_array(3),status)
  enddo

! create extension QC
! contains reliability flags

! break if there were errors

  if (status.ne.0) then
    call ftgerr(status,cfitsioerror)
    print *,gettime(),"FITS output error: ",trim(filename),": ",status,cfitsioerror
    call exit(1)
  endif

! close

  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine write_fits

subroutine create_output_arrays(all_results,runs,resultprocessingarray,resultprocessingtext)
!define arrays with links to all the data that needs processing.
!extinction, diagnostics, cel abundances, orl abundances, strong line
!abundances, adfs

  implicit none
  real(kind=dp), dimension(:,:), allocatable :: resultprocessingarray
  character(len=40), dimension(:,:), allocatable :: resultprocessingtext
  character(len=35) :: extinction_format, diagnostic_format, diagnostic_ratio_format, abundances_format, adf_format
  type(resultarray), dimension(:) :: all_results
  integer :: runs

  extinction_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
  diagnostic_format = "(X,A,F9.1,SP,F9.1,F9.1,S)"
  diagnostic_ratio_format = "(X,A,F9.3,SP,F9.3,F9.3,S)"
  abundances_format = "(X,A,ES14.3,SP,ES14.3,ES14.3,S)"
  adf_format = "(X,A,F8.2,SP,F8.2,F8.2,S)"

  allocate(resultprocessingarray(164,runs))
  resultprocessingarray=0d0
  allocate(resultprocessingtext(164,4))

!extinction

  resultprocessingarray(1,:) = all_results%cHb_ha
  resultprocessingtext(1,:) = (/"c(Hb) (Ha/Hb)                      ","c(H$\beta)$ (H$\alpha/H$\beta)     ", extinction_format, "chb_ha_hb                          "/)
  resultprocessingarray(2,:) = all_results%cHb_hg
  resultprocessingtext(2,:) = (/"c(Hb) (Hg/Hb)                      ","c(H$\beta)$ (H$\gamma/H$\beta)     ", extinction_format, "chb_hg_hb                          "/)
  resultprocessingarray(3,:) = all_results%cHb_hd
  resultprocessingtext(3,:) = (/"c(Hb) (Hd/Hb)                      ","c(H$\beta)$ (H$\delta/H$\beta)     ", extinction_format, "chb_hd_hb                          "/)
  resultprocessingarray(4,:) = all_results%mean_cHb
  resultprocessingtext(4,:) = (/"mean c(Hb)                         ","c(H$\beta)$                        ", extinction_format, "mean_chb                           "/)

!diagnostics
!low density

  resultprocessingarray(5,:) = all_results%oii_density
  resultprocessingtext(5,:) = (/"[OII] density                      ","{}[O~{\sc ii}] density             ", diagnostic_format, "density_oii                        "/)
  resultprocessingarray(6,:) = all_results%oii_density_ratio
  resultprocessingtext(6,:) = (/"[OII] 3729/3726 ratio              ","{}[O~{\sc ii}] 3729/3726 ratio     ", diagnostic_ratio_format, "density_oii_ratio                  "/)
  resultprocessingarray(7,:) = all_results%SII_density
  resultprocessingtext(7,:) = (/"[SII] density                      ","{}[S~{\sc ii}] density             ", diagnostic_format, "density_sii                        "/)
  resultprocessingarray(8,:) = all_results%sii_density_ratio
  resultprocessingtext(8,:) = (/"[SII] 6731/6717 ratio              ","{}[S~{\sc ii}] 6731/6717 ratio     ", diagnostic_ratio_format, "density_sii_ratio                  "/)
  resultprocessingarray(9,:) = all_results%low_density
  resultprocessingtext(9,:) = (/"Low ionisation density             ","Low ionisation density             ", diagnostic_format, "density_low                        "/)

!low temperature

  resultprocessingarray(10,:) = all_results%oii_temp
  resultprocessingtext(10,:) = (/"[OII] temperature                  ","{}[O~{\sc ii}] temperature         ", diagnostic_format, "temp_oii                           "/)
  resultprocessingarray(11,:) = all_results%oii_temp_ratio
  resultprocessingtext(11,:) = (/"[OII] (7320+7330)/(3726+3729) ratio","{}[O~{\sc ii}] 7320+7330/3726+3729 ", diagnostic_ratio_format, "temp_oii_ratio                     "/)
  resultprocessingarray(12,:) = all_results%SII_temp
  resultprocessingtext(12,:) = (/"[SII] temperature                  ","{}[S~{\sc ii}] temperature         ", diagnostic_format, "temp_sii                           "/)
  resultprocessingarray(13,:) = all_results%sii_temp_ratio
  resultprocessingtext(13,:) = (/"[SII] (6717+6731)/(4068+4076) ratio","{}[S~{\sc ii}] 6717+6731/4068+4076 ", diagnostic_ratio_format, "temp_sii_ratio                     "/)
  resultprocessingarray(14,:) = all_results%NII_temp
  resultprocessingtext(14,:) = (/"[NII] temperature                  ","{}[N~{\sc ii}] temperature         ", diagnostic_format, "temp_nii                           "/)
  resultprocessingarray(15,:) = all_results%nii_temp_ratio
  resultprocessingtext(15,:) = (/"[NII] (6548+6584)/5754 ratio       ","{}[N~{\sc ii}] 6548+6584/5754 ratio", diagnostic_ratio_format, "temp_nii_ratio                     "/)
  resultprocessingarray(16,:) = all_results%OI_temp
  resultprocessingtext(16,:) = (/"[OI] temperature                   ","{}[O~{\sc i}] temperature          ", diagnostic_format, "temp_oi                            "/)
  resultprocessingarray(17,:) = all_results%oi_temp_ratio
  resultprocessingtext(17,:) = (/"[OI] (6300+6363)/5577 ratio        ","{}[O~{\sc i}] 6300+6363/5577 ratio ", diagnostic_ratio_format, "temp_oi_ratio                      "/)
  resultprocessingarray(18,:) = all_results%CI_temp
  resultprocessingtext(18,:) = (/"[CI] temperature                   ","{}[C~{\sc i}] temperature          ", diagnostic_format, "temp_ci                            "/)
  resultprocessingarray(19,:) = all_results%ci_temp_ratio
  resultprocessingtext(19,:) = (/"[CI] (9824+9850)/8727 ratio        ","{}[C~{\sc i}] 9824+9850/8727 ratio ", diagnostic_ratio_format, "temp_ci_ratio                      "/)
  resultprocessingarray(20,:) = all_results%low_temp
  resultprocessingtext(20,:) = (/"Low ionisation temperature         ","Low ionisation temperature         ", diagnostic_format, "temp_low                           "/)


!medium density

  resultprocessingarray(21,:) = all_results%cliii_density
  resultprocessingtext(21,:) = (/"[ClIII] density                    ","{}[Cl~{\sc iii}] density           ", diagnostic_format, "density_cliii                      "/)
  resultprocessingarray(22,:) = all_results%cliii_density_ratio
  resultprocessingtext(22,:) = (/"[ClIII] 5537/5517 ratio            ","{}[Cl~{\sc iii}] 5537/5517 ratio   ", diagnostic_ratio_format, "density_cliii_ratio                "/)
  resultprocessingarray(23,:) = all_results%ArIV_density
  resultprocessingtext(23,:) = (/"[ArIV] density                     ","{}[Ar~{\sc iv}] density            ", diagnostic_format, "density_ariv                       "/)
  resultprocessingarray(24,:) = all_results%ariv_density_ratio
  resultprocessingtext(24,:) = (/"[ArIV] 4740/4711 ratio             ","{}[Ar~{\sc iv}] 4740/4711 ratio    ", diagnostic_ratio_format, "density_sii_ratio                  "/)
  resultprocessingarray(25,:) = all_results%CIII_density
  resultprocessingtext(25,:) = (/"[CIII] density                     ","{}[C~{\sc iii}] density            ", diagnostic_format, "density_ciii                       "/)
  resultprocessingarray(26,:) = all_results%ciii_density_ratio
  resultprocessingtext(26,:) = (/"[CIII] 1909/1907 ratio             ","{}[C~{\sc iii}] 1909/1907 ratio    ", diagnostic_ratio_format, "density_ciii_ratio                 "/)
  resultprocessingarray(27,:) = all_results%OIII_IR_density
  resultprocessingtext(27,:) = (/"[OIII] IR density                  ","{}[O~{\sc iii}] IR density         ", diagnostic_format, "density_oiii_ir                    "/)
  resultprocessingarray(28,:) = all_results%oiii_ir_density_ratio
  resultprocessingtext(28,:) = (/"[OIII] IR 88um/52um ratio          ","{}[O~{\sc iii}] IR 88um/52um ratio ", diagnostic_ratio_format, "density_oiii_ir_ratio              "/)
  resultprocessingarray(29,:) = all_results%SIII_IR_density
  resultprocessingtext(29,:) = (/"[SIII] IR density                  ","{}[S~{\sc iii}] IR density         ", diagnostic_format, "density_siii_ir                    "/)
  resultprocessingarray(30,:) = all_results%siii_ir_density_ratio
  resultprocessingtext(30,:) = (/"[SIII] IR 33um/18um ratio          ","{}[S~{\sc iii}] IR 33um/18um ratio ", diagnostic_ratio_format, "density_siii_ir_ratio              "/)
  resultprocessingarray(31,:) = all_results%ArIII_IR_density
  resultprocessingtext(31,:) = (/"[ArIII] IR density                 ","{}[Ar~{\sc iii}] IR density        ", diagnostic_format, "density_ariii_ir                   "/)
  resultprocessingarray(32,:) = all_results%ariii_ir_density_ratio
  resultprocessingtext(32,:) = (/"[ArIII] 9um/22um ratio             ","{}[Ar~{\sc iii}] IR 9um/22um ratio ", diagnostic_ratio_format, "density_ariii_ir_ratio             "/)
  resultprocessingarray(33,:) = all_results%NeIII_IR_density
  resultprocessingtext(33,:) = (/"[NeIII] IR density                 ","{}[Ne~{\sc iii}] IR density        ", diagnostic_format, "density_neiii_ir                   "/)
  resultprocessingarray(34,:) = all_results%neiii_ir_density_ratio
  resultprocessingtext(34,:) = (/"[Ne III] 15um/36um ratio           ","{}[Ne~{\sc iii}] IR 15um/36um ratio", diagnostic_ratio_format, "density_neiii_ir_ratio             "/)
  resultprocessingarray(35,:) = all_results%med_density
  resultprocessingtext(35,:) = (/"Medium ionisation density          ","Medium ionisation density          ", diagnostic_format, "density_med                        "/)

! medium temperature

  resultprocessingarray(36,:) = all_results%OIII_temp
  resultprocessingtext(36,:) = (/"[OIII] temperature                 ","{}[O~{\sc iii}] temperature        ", diagnostic_format, "temp_oiii                          "/)
  resultprocessingarray(37,:) = all_results%oiii_temp_ratio
  resultprocessingtext(37,:) = (/"[OIII] (4959+5007)/4363 ratio      ","{}[O~{\sc iii}] 4959+5007/4363     ", diagnostic_ratio_format, "temp_oiii_ratio                    "/)
  resultprocessingarray(38,:) = all_results%OIII_IR_temp
  resultprocessingtext(38,:) = (/"[OIII] IR temperature              ","{}[O~{\sc iii}] IR temperature     ", diagnostic_format, "temp_oiii_ir                       "/)
  resultprocessingarray(39,:) = all_results%oiii_ir_temp_ratio
  resultprocessingtext(39,:) = (/"[OIII] (4959+5007)/52um ratio      ","{}[O~{\sc iii}] 4959+5007/52um     ", diagnostic_ratio_format, "temp_oiii_ir_ratio                 "/)
  resultprocessingarray(40,:) = all_results%OIII_UV_temp
  resultprocessingtext(40,:) = (/"[OIII] UV temperature              ","{}[O~{\sc iii}] UV temperature     ", diagnostic_format, "temp_oiii_uv                       "/)
  resultprocessingarray(41,:) = all_results%oiii_uv_temp_ratio
  resultprocessingtext(41,:) = (/"[OIII] (4959+5007)/1666 ratio      ","{}[O~{\sc iii}] 4959+5007/1666     ", diagnostic_ratio_format, "temp_oiii_uv_ratio                 "/)
  resultprocessingarray(42,:) = all_results%NeIII_temp
  resultprocessingtext(42,:) = (/"[NeIII] temperature                ","{}[Ne~{\sc iii}] temperature       ", diagnostic_format, "temp_neiii                         "/)
  resultprocessingarray(43,:) = all_results%neiii_temp_ratio
  resultprocessingtext(43,:) = (/"[NeIII] (3868+3967)/3342 ratio     ","{}[Ne~{\sc iii}] 3868+3967/3342    ", diagnostic_ratio_format, "temp_neiii_ratio                   "/)
  resultprocessingarray(44,:) = all_results%NeIII_IR_temp
  resultprocessingtext(44,:) = (/"[NeIII] IR temperature             ","{}[Ne~{\sc iii}] IR temperature    ", diagnostic_format, "temp_neiii_ir                      "/)
  resultprocessingarray(45,:) = all_results%neiii_ir_temp_ratio
  resultprocessingtext(45,:) = (/"[NeIII] (3868+3967)/15um ratio     ","{}[Ne~{\sc iii}] 3868+3967/15um    ", diagnostic_ratio_format, "temp_neiii_ir_ratio                "/)
  resultprocessingarray(46,:) = all_results%ArIII_temp
  resultprocessingtext(46,:) = (/"[ArIII] temperature                ","{}[Ar~{\sc iii}] temperature       ", diagnostic_format, "temp_ariii                         "/)
  resultprocessingarray(47,:) = all_results%ariii_temp_ratio
  resultprocessingtext(47,:) = (/"[ArIII] (7135+7751)/5192 ratio     ","{}[Ar~{\sc iii}] 7135+7751/5192    ", diagnostic_ratio_format, "temp_ariii_ratio                   "/)
  resultprocessingarray(48,:) = all_results%SIII_temp
  resultprocessingtext(48,:) = (/"[SIII] temperature                 ","{}[S~{\sc iii}] temperature        ", diagnostic_format, "temp_siii                          "/)
  resultprocessingarray(49,:) = all_results%siii_temp_ratio
  resultprocessingtext(49,:) = (/"[SIII] (9069+9531)/6312 ratio      ","{}[S~{\sc iii}] 9069+9531/6312     ", diagnostic_ratio_format, "temp_siii_ratio                    "/)
  resultprocessingarray(50,:) = all_results%med_temp
  resultprocessingtext(50,:) = (/"Medium ionisation temperature      ","Medium ionisation temperature      ", diagnostic_format, "temp_med                           "/)

!high density

  resultprocessingarray(51,:) = all_results%neiv_density
  resultprocessingtext(51,:) = (/"[NeIV] density                     ","{}[Ne~{\sc iv}] density            ", diagnostic_format, "density_neiv                       "/)
  resultprocessingarray(52,:) = all_results%neiv_density_ratio
  resultprocessingtext(52,:) = (/"[NeIV] 2425/2423 ratio             ","{}[Ne~{\sc iv}] 2425/2423 ratio    ", diagnostic_ratio_format, "density_neiv_ratio                 "/)
  resultprocessingarray(53,:) = all_results%high_density
  resultprocessingtext(53,:) = (/"High ionisation density            ","High ionisation density            ", diagnostic_format, "density_high                       "/)

!high temperature

  resultprocessingarray(54,:) = all_results%ArV_temp
  resultprocessingtext(54,:) = (/"[ArV] temperature                  ","{}[Ar~{\sc v}] temperature         ", diagnostic_format, "temp_arv                           "/)
  resultprocessingarray(55,:) = all_results%arv_temp_ratio
  resultprocessingtext(55,:) = (/"[ArV] (6435+7005)/4625 ratio       ","{}[Ar~{\sc v}] 6435+7005/4625 ratio", diagnostic_ratio_format, "temp_arv_ratio                     "/)
  resultprocessingarray(56,:) = all_results%NeV_temp
  resultprocessingtext(56,:) = (/"[NeV] temperature                  ","{}[Ne~{\sc v}] temperature         ", diagnostic_format, "temp_nev                           "/)
  resultprocessingarray(57,:) = all_results%nev_temp_ratio
  resultprocessingtext(57,:) = (/"[NeV] (3426+3345)/2975 ratio       ","{}[Ne~{\sc v}] 3426+3345/2975 ratio", diagnostic_ratio_format, "temp_nev_ratio                     "/)
  resultprocessingarray(58,:) = all_results%high_temp
  resultprocessingtext(58,:) = (/"High temperature                   ","High temperature                   ", diagnostic_format, "temp_high                          "/)

!Recombination line diagnostics

  resultprocessingarray(59,:) = all_results%Balmer_jump_temp
  resultprocessingtext(59,:) = (/"Balmer jump temperature            ","Balmer jump temperature            ", diagnostic_format, "temp_balmerjump                    "/)

  resultprocessingarray(60,:) = all_results%Paschen_jump_temp
  resultprocessingtext(60,:) = (/"Paschen jump temperature           ","Paschen jump temperature           ", diagnostic_format, "temp_paschenjump                   "/)

  resultprocessingarray(61,:) = all_results%Balmerdec_density
  resultprocessingtext(61,:) = (/"Balmer decrement density           ","Balmer decrement density           ", diagnostic_format, "density_balmerdec                  "/)

  resultprocessingarray(62,:) = all_results%Paschendec_density
  resultprocessingtext(62,:) = (/"Paschen decrement density          ","Paschen decrement density          ", diagnostic_format, "density_paschendec                 "/)

  resultprocessingarray(63,:) = all_results%te_5876_4471
  resultprocessingtext(63,:) = (/"He I temperature (5876/4471)       ","He I temperature (5876/4471)       ", diagnostic_format, "temp_he_5876_4471                  "/)

  resultprocessingarray(64,:) = all_results%ratio_5876_4471
  resultprocessingtext(64,:) = (/"He I 5876/4471 ratio               ","He I 5876/4471 ratio               ", diagnostic_ratio_format, "temp_he_5876_4471_ratio            "/)

  resultprocessingarray(65,:) = all_results%te_6678_4471
  resultprocessingtext(65,:) = (/"He I temperature (6678/4471)       ","He I temperature (6678/4471)       ", diagnostic_format, "temp_he_6678_4471                  "/)

  resultprocessingarray(66,:) = all_results%ratio_6678_4471
  resultprocessingtext(66,:) = (/"He I 6678/4471 ratio               ","He I 6678/4471 ratio               ", diagnostic_ratio_format, "temp_he_6678_4471_ratio            "/)

  resultprocessingarray(67,:) = all_results%oii_rl_R1
  resultprocessingtext(67,:) = (/"OII 4649/4089 ratio                ","OII 4649/4089 ratio                ", diagnostic_ratio_format, "temp_oiirls_R1                     "/)

  resultprocessingarray(68,:) = all_results%oii_rl_R2
  resultprocessingtext(68,:) = (/"OII 4649/4662 ratio                ","OII 4649/4662 ratio                ", diagnostic_ratio_format, "temp_oiirls_R2                     "/)

  resultprocessingarray(69,:) = all_results%oii_te
  resultprocessingtext(69,:) = (/"OII temperature                    ","OII temperature                    ", diagnostic_format, "temp_oiirls                        "/)

  resultprocessingarray(70,:) = all_results%oii_ne
  resultprocessingtext(70,:) = (/"OII density                        ","OII density                        ", diagnostic_format, "density_oiirls                     "/)

! recombination contribution to CELs

  resultprocessingarray(71,:) = all_results%nii5754recCEL
  resultprocessingtext(71,:) = (/"NII 5754 assuming N2+ from CELs    ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_CELabund        "/)
  resultprocessingarray(72,:) = all_results%nii5754recRL
  resultprocessingtext(72,:) = (/"NII 5754 assuming N2+ from RLs     ","N~{\sc ii}5754$_R$ (N$^{2+}$ CELs) ", adf_format, "recombination_5754_RLabund         "/)
  resultprocessingarray(73,:) = all_results%oii7325recCEL
  resultprocessingtext(73,:) = (/"OII 7320,30 assuming O2+ from CELs ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_CELabund        "/)
  resultprocessingarray(74,:) = all_results%oii7325recRL
  resultprocessingtext(74,:) = (/"OII 7320,30 assuming O2+ from RLs  ","O~{\sc ii}7325$_R$ (O$^{2+}$ CELs) ", adf_format, "recombination_7325_RLabund         "/)
  resultprocessingarray(75,:) = all_results%oiii4363recCEL
  resultprocessingtext(75,:) = (/"OIII 4363 assuming O3+ from CELs   ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_CELabund        "/)
  resultprocessingarray(76,:) = all_results%oiii4363recRL
  resultprocessingtext(76,:) = (/"OIII 4363 assuming O3+ from RLs    ","O~{\sc ii}4363$_R$ (O$^{3+}$ CELs) ", adf_format, "recombination_4363_RLabund         "/)

!CEL abundances

  resultprocessingarray(77,:) = all_results%NC_abund_CEL
  resultprocessingtext(77,:) = (/"C0/H                               ","C$^{0}$/H                          ", abundances_format, "abund_cel_nc                       "/)
  resultprocessingarray(78,:) = all_results%cii_abund_CEL
  resultprocessingtext(78,:) = (/"C+/H                               ","C$^{+}$/H                          ", abundances_format, "abund_cel_cii                      "/)
  resultprocessingarray(79,:) = all_results%ciii_abund_CEL
  resultprocessingtext(79,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_cel_ciii                     "/)
  resultprocessingarray(80,:) = all_results%civ_abund_CEL
  resultprocessingtext(80,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_cel_civ                      "/)
  resultprocessingarray(81,:) = all_results%c_icf_CEL
  resultprocessingtext(81,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_cel_c                          "/)
  resultprocessingarray(82,:) = all_results%C_abund_CEL
  resultprocessingtext(82,:) = (/"C/H                                ","C$^{}$/H                           ", abundances_format, "abund_cel_c                        "/)
  resultprocessingarray(83,:) = all_results%nii_abund_CEL
  resultprocessingtext(83,:) = (/"N+/H                               ","N$^{+}$/H                          ", abundances_format, "abund_cel_nii                      "/)
  resultprocessingarray(84,:) = all_results%niii_abund_CEL
  resultprocessingtext(84,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_cel_niii                     "/)
  resultprocessingarray(85,:) = all_results%niv_abund_CEL
  resultprocessingtext(85,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_cel_niv                      "/)
  resultprocessingarray(86,:) = all_results%nv_abund_CEL
  resultprocessingtext(86,:) = (/"N4+/H                              ","N$^{4+}$/H                         ", abundances_format, "abund_cel_nv                       "/)
  resultprocessingarray(87,:) = all_results%n_icf_CEL
  resultprocessingtext(87,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_cel_n                          "/)
  resultprocessingarray(88,:) = all_results%N_abund_CEL
  resultprocessingtext(88,:) = (/"N/H                                ","N$^{}$/H                           ", abundances_format, "abund_cel_n                        "/)
  resultprocessingarray(89,:) = all_results%NO_abund_CEL
  resultprocessingtext(89,:) = (/"O0/H                               ","O$^{0}$/H                          ", abundances_format, "abund_cel_no                       "/)
  resultprocessingarray(90,:) = all_results%Oii_abund_CEL
  resultprocessingtext(90,:) = (/"O+/H                               ","O$^{+}$/H                          ", abundances_format, "abund_cel_oii                      "/)
  resultprocessingarray(91,:) = all_results%Oiii_abund_CEL
  resultprocessingtext(91,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_cel_oiii                     "/)
  resultprocessingarray(92,:) = all_results%Oiv_abund_CEL
  resultprocessingtext(92,:) = (/"O3+/H                              ","O$^{3+}$/H                         ", abundances_format, "abund_cel_oiv                      "/)
  resultprocessingarray(93,:) = all_results%o_icf_CEL
  resultprocessingtext(93,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_cel_o                          "/)
  resultprocessingarray(94,:) = all_results%O_abund_CEL
  resultprocessingtext(94,:) = (/"O/H                                ","O$^{}$/H                           ", abundances_format, "abund_cel_o                        "/)
  resultprocessingarray(95,:) = all_results%NeII_abund_CEL
  resultprocessingtext(95,:) = (/"Ne+/H                              ","Ne$^{+}$/H                         ", abundances_format, "abund_cel_neii                     "/)
  resultprocessingarray(96,:) = all_results%NeIII_abund_CEL
  resultprocessingtext(96,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_cel_neiii                    "/)
  resultprocessingarray(97,:) = all_results%NeIV_abund_CEL
  resultprocessingtext(97,:) = (/"Ne3+/H                             ","Ne$^{3+}$/H                        ", abundances_format, "abund_cel_neiv                     "/)
  resultprocessingarray(98,:) = all_results%NeV_abund_CEL
  resultprocessingtext(98,:) = (/"Ne4+/H                             ","Ne$^{4+}$/H                        ", abundances_format, "abund_cel_nev                      "/)
  resultprocessingarray(99,:) = all_results%ne_icf_CEL
  resultprocessingtext(99,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_cel_ne                         "/)
  resultprocessingarray(100,:) = all_results%Ne_abund_CEL
  resultprocessingtext(100,:) = (/"Ne/H                               ","Ne$^{}$/H                          ", abundances_format, "abund_cel_ne                       "/)
  resultprocessingarray(101,:) = all_results%ArIII_abund_CEL
  resultprocessingtext(101,:) = (/"Ar2+/H                             ","Ar$^{2+}$/H                        ", abundances_format, "abund_cel_ariii                    "/)
  resultprocessingarray(102,:) = all_results%ArIV_abund_CEL
  resultprocessingtext(102,:) = (/"Ar3+/H                             ","Ar$^{3+}$/H                        ", abundances_format, "abund_cel_ariv                     "/)
  resultprocessingarray(103,:) = all_results%ArV_abund_CEL
  resultprocessingtext(103,:) = (/"Ar4+/H                             ","Ar$^{4+}$/H                        ", abundances_format, "abund_cel_arv                      "/)
  resultprocessingarray(104,:) = all_results%ar_icf_CEL
  resultprocessingtext(104,:) = (/"icf(Ar)                            ","icf(Ar)                            ", abundances_format, "icf_cel_ar                         "/)
  resultprocessingarray(105,:) = all_results%Ar_abund_CEL
  resultprocessingtext(105,:) = (/"Ar/H                               ","Ar$^{}$/H                          ", abundances_format, "abund_cel_ar                       "/)
  resultprocessingarray(106,:) = all_results%SII_abund_CEL
  resultprocessingtext(106,:) = (/"S+/H                               ","S$^{+}$/H                          ", abundances_format, "abund_cel_sii                      "/)
  resultprocessingarray(107,:) = all_results%SIII_abund_CEL
  resultprocessingtext(107,:) = (/"S2+/H                              ","S$^{2+}$/H                         ", abundances_format, "abund_cel_siii                     "/)
  resultprocessingarray(108,:) = all_results%s_icf_CEL
  resultprocessingtext(108,:) = (/"icf(S)                             ","icf(S)                             ", abundances_format, "icf_cel_s                          "/)
  resultprocessingarray(109,:) = all_results%S_abund_CEL
  resultprocessingtext(109,:) = (/"S/H                                ","S$^{}$/H                           ", abundances_format, "abund_cel_s                        "/)

  resultprocessingarray(110,:) = all_results%ClII_abund_CEL
  resultprocessingtext(110,:) = (/"Cl+/H                              ","Cl$^{+}$/H                         ", abundances_format, "abund_cel_clii                     "/)
  resultprocessingarray(111,:) = all_results%ClIII_abund_CEL
  resultprocessingtext(111,:) = (/"Cl2+/H                             ","Cl$^{2+}$/H                        ", abundances_format, "abund_cel_cliii                    "/)
  resultprocessingarray(112,:) = all_results%ClIV_abund_CEL
  resultprocessingtext(112,:) = (/"Cl3+/H                             ","Cl$^{3+}$/H                        ", abundances_format, "abund_cel_cliv                     "/)

  resultprocessingarray(113,:) = all_results%cl_icf_CEL
  resultprocessingtext(113,:) = (/"icf(Cl)                            ","icf(Cl)                            ", abundances_format, "icf_cel_cl                         "/)
  resultprocessingarray(114,:) = all_results%Cl_abund_CEL
  resultprocessingtext(114,:) = (/"Cl/H                               ","Cl$^{}$/H                          ", abundances_format, "abund_cel_cl                       "/)

!ORL abundances

  resultprocessingarray(115,:) = all_results%Hei_abund_ORL
  resultprocessingtext(115,:) = (/"He+/H                              ","He$^{+}$/H                         ", abundances_format, "abund_orl_hei                      "/)
  resultprocessingarray(116,:) = all_results%Heii_abund_ORL
  resultprocessingtext(116,:) = (/"He2+/H                             ","He$^{2+}$/H                        ", abundances_format, "abund_orl_heii                     "/)
  resultprocessingarray(117,:) = all_results%He_abund_ORL
  resultprocessingtext(117,:) = (/"He/H                               ","He/H                               ", abundances_format, "abund_orl_he                       "/)
  resultprocessingarray(118,:) = all_results%Cii_abund_ORL
  resultprocessingtext(118,:) = (/"C2+/H                              ","C$^{2+}$/H                         ", abundances_format, "abund_orl_cii                      "/)
  resultprocessingarray(119,:) = all_results%Ciii_abund_ORL
  resultprocessingtext(119,:) = (/"C3+/H                              ","C$^{3+}$/H                         ", abundances_format, "abund_orl_ciii                     "/)
  resultprocessingarray(120,:) = all_results%c_icf_ORL
  resultprocessingtext(120,:) = (/"icf(C)                             ","icf(C)                             ", abundances_format, "icf_orl_c                          "/)
  resultprocessingarray(121,:) = all_results%C_abund_ORL
  resultprocessingtext(121,:) = (/"C/H                                ","C/H                                ", abundances_format, "abund_orl_c                        "/)
  resultprocessingarray(122,:) = all_results%Nii_v3_abund_ORL
  resultprocessingtext(122,:) = (/"N2+/H (V3)                         ","N$^{2+}$/H (V3)                    ", abundances_format, "abund_orl_nii_v3                   "/)
  resultprocessingarray(123,:) = all_results%Nii_v5_abund_ORL
  resultprocessingtext(123,:) = (/"N2+/H (V5)                         ","N$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_nii_v5                   "/)
  resultprocessingarray(124,:) = all_results%Nii_v8_abund_ORL
  resultprocessingtext(124,:) = (/"N2+/H (V8)                         ","N$^{2+}$/H (V8)                    ", abundances_format, "abund_orl_nii_v8                   "/)
  resultprocessingarray(125,:) = all_results%Nii_v12_abund_ORL
  resultprocessingtext(125,:) = (/"N2+/H (V12)                        ","N$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_nii_v12                  "/)
  resultprocessingarray(126,:) = all_results%Nii_v20_abund_ORL
  resultprocessingtext(126,:) = (/"N2+/H (V20)                        ","N$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_nii_v20                  "/)
  resultprocessingarray(127,:) = all_results%Nii_v28_abund_ORL
  resultprocessingtext(127,:) = (/"N2+/H (V28)                        ","N$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_nii_v28                  "/)
  resultprocessingarray(128,:) = all_results%Nii_3d4f_abund_ORL
  resultprocessingtext(128,:) = (/"N2+/H (3d-4f)                      ","N$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_nii_3d4f                 "/)
  resultprocessingarray(129,:) = all_results%Nii_abund_ORL
  resultprocessingtext(129,:) = (/"N2+/H                              ","N$^{2+}$/H                         ", abundances_format, "abund_orl_nii                      "/)
  resultprocessingarray(130,:) = all_results%Niii_abund_ORL
  resultprocessingtext(130,:) = (/"N3+/H                              ","N$^{3+}$/H                         ", abundances_format, "abund_orl_niii                     "/)
  resultprocessingarray(131,:) = all_results%n_icf_ORL
  resultprocessingtext(131,:) = (/"icf(N)                             ","icf(N)                             ", abundances_format, "icf_orl_n                          "/)
  resultprocessingarray(132,:) = all_results%N_abund_ORL
  resultprocessingtext(132,:) = (/"N/H                                ","N/H                                ", abundances_format, "abund_orl_n                        "/)
  resultprocessingarray(133,:) = all_results%Oii_v1_abund_ORL
  resultprocessingtext(133,:) = (/"O2+/H (V1)                         ","O$^{2+}$/H (V1)                    ", abundances_format, "abund_orl_oii_v1                   "/)
  resultprocessingarray(134,:) = all_results%Oii_v2_abund_ORL
  resultprocessingtext(134,:) = (/"O2+/H (V2)                         ","O$^{2+}$/H (V2)                    ", abundances_format, "abund_orl_oii_v2                   "/)
  resultprocessingarray(135,:) = all_results%Oii_v5_abund_ORL
  resultprocessingtext(135,:) = (/"O2+/H (V5)                         ","O$^{2+}$/H (V5)                    ", abundances_format, "abund_orl_oii_v5                   "/)
  resultprocessingarray(136,:) = all_results%Oii_v10_abund_ORL
  resultprocessingtext(136,:) = (/"O2+/H (V10)                        ","O$^{2+}$/H (V10)                   ", abundances_format, "abund_orl_oii_v10                  "/)
  resultprocessingarray(137,:) = all_results%Oii_v11_abund_ORL
  resultprocessingtext(137,:) = (/"O2+/H (V11)                        ","O$^{2+}$/H (V11)                   ", abundances_format, "abund_orl_oii_v11                  "/)
  resultprocessingarray(138,:) = all_results%Oii_v12_abund_ORL
  resultprocessingtext(138,:) = (/"O2+/H (V12)                        ","O$^{2+}$/H (V12)                   ", abundances_format, "abund_orl_oii_v12                  "/)
  resultprocessingarray(139,:) = all_results%Oii_v19_abund_ORL
  resultprocessingtext(139,:) = (/"O2+/H (V19)                        ","O$^{2+}$/H (V19)                   ", abundances_format, "abund_orl_oii_v19                  "/)
  resultprocessingarray(140,:) = all_results%Oii_v20_abund_ORL
  resultprocessingtext(140,:) = (/"O2+/H (V20)                        ","O$^{2+}$/H (V20)                   ", abundances_format, "abund_orl_oii_v20                  "/)
  resultprocessingarray(141,:) = all_results%Oii_v25_abund_ORL
  resultprocessingtext(141,:) = (/"O2+/H (V25)                        ","O$^{2+}$/H (V25)                   ", abundances_format, "abund_orl_oii_v25                  "/)
  resultprocessingarray(142,:) = all_results%Oii_v28_abund_ORL
  resultprocessingtext(142,:) = (/"O2+/H (V28)                        ","O$^{2+}$/H (V28)                   ", abundances_format, "abund_orl_oii_v28                  "/)
  resultprocessingarray(143,:) = all_results%Oii_v33_abund_ORL
  resultprocessingtext(143,:) = (/"O2+/H (V33)                        ","O$^{2+}$/H (V33)                   ", abundances_format, "abund_orl_oii_v33                  "/)
  resultprocessingarray(144,:) = all_results%Oii_3d4f_abund_ORL
  resultprocessingtext(144,:) = (/"O2+/H (3d-4f)                      ","O$^{2+}$/H (3d-4f)                 ", abundances_format, "abund_orl_oii_3d4f                 "/)
  resultprocessingarray(145,:) = all_results%Oii_abund_ORL
  resultprocessingtext(145,:) = (/"O2+/H                              ","O$^{2+}$/H                         ", abundances_format, "abund_orl_oii                      "/)
  resultprocessingarray(146,:) = all_results%O_icf_ORL
  resultprocessingtext(146,:) = (/"icf(O)                             ","icf(O)                             ", abundances_format, "icf_orl_o                          "/)
  resultprocessingarray(147,:) = all_results%O_abund_ORL
  resultprocessingtext(147,:) = (/"O/H                                ","O/H                                ", abundances_format, "abund_orl_o                        "/)
  resultprocessingarray(148,:) = all_results%Neii_abund_ORL
  resultprocessingtext(148,:) = (/"Ne2+/H                             ","Ne$^{2+}$/H                        ", abundances_format, "abund_orl_neii                     "/)
  resultprocessingarray(149,:) = all_results%ne_icf_ORL
  resultprocessingtext(149,:) = (/"icf(Ne)                            ","icf(Ne)                            ", abundances_format, "icf_orl_ne                         "/)
  resultprocessingarray(150,:) = all_results%Ne_abund_ORL
  resultprocessingtext(150,:) = (/"Ne/H                               ","Ne/H                               ", abundances_format, "abund_orl_ne                       "/)

!strong line abundances

  resultprocessingarray(151,:) = all_results%O_R23_upper
  resultprocessingtext(151,:) = (/"O/H (R23 upper)                    ","O/H (R23 upper)                    ", abundances_format, "abund_o_r23_upper                  "/)
  resultprocessingarray(152,:) = all_results%O_R23_lower
  resultprocessingtext(152,:) = (/"O/H (R23 lower)                    ","O/H (R23 lower)                    ", abundances_format, "abund_o_r23_lower                  "/)
  resultprocessingarray(153,:) = all_results%O_N2
  resultprocessingtext(153,:) = (/"O/H (N2)                           ","O/H (N2)                           ", abundances_format, "abund_o_n2                         "/)
  resultprocessingarray(154,:) = all_results%O_O3N2
  resultprocessingtext(154,:) = (/"O/H (O3N2)                         ","O/H (O3N2)                         ", abundances_format, "abund_o_o3n2                       "/)
  resultprocessingarray(155,:) = all_results%O_Ar3O3
  resultprocessingtext(155,:) = (/"O/H (Ar3O3)                        ","O/H (Ar3O3)                        ", abundances_format, "abund_o_ar3o3                      "/)
  resultprocessingarray(156,:) = all_results%O_S3O3
  resultprocessingtext(156,:) = (/"O/H (S3O3)                         ","O/H (S3O3)                         ", abundances_format, "abund_o_s3o3                       "/)

!adfs

  resultprocessingarray(157,:) = all_results%adf_o2plus
  resultprocessingtext(157,:) = (/"adf (O2+/H)                        ","adf (O$^{2+}$/H)                   ", adf_format, "adf_o2plus                         "/)
  resultprocessingarray(158,:) = all_results%adf_o
  resultprocessingtext(158,:) = (/"adf (O/H)                          ","adf (O/H)                          ", adf_format, "adf_o                              "/)
  resultprocessingarray(159,:) = all_results%adf_n2plus
  resultprocessingtext(159,:) = (/"adf (N2+/H)                        ","adf (N$^{2+}$/H)                   ", adf_format, "adf_n2plus                         "/)
  resultprocessingarray(160,:) = all_results%adf_n
  resultprocessingtext(160,:) = (/"adf (N/H)                          ","adf (N/H)                          ", adf_format, "adf_n                              "/)
  resultprocessingarray(161,:) = all_results%adf_c2plus
  resultprocessingtext(161,:) = (/"adf (C2+/H)                        ","adf (C$^{2+}$/H)                   ", adf_format, "adf_c2plus                         "/)
  resultprocessingarray(162,:) = all_results%adf_c
  resultprocessingtext(162,:) = (/"adf (C/H)                          ","adf (C/H)                          ", adf_format, "adf_c                              "/)
  resultprocessingarray(163,:) = all_results%adf_ne2plus
  resultprocessingtext(163,:) = (/"adf (Ne2+/H)                       ","adf (Ne$^{2+}$/H)                  ", adf_format, "adf_ne2plus                        "/)
  resultprocessingarray(164,:) = all_results%adf_ne
  resultprocessingtext(164,:) = (/"adf (Ne/H)                         ","adf (Ne/H)                         ", adf_format, "adf_ne                             "/)

end subroutine create_output_arrays

end module mod_output
