!functions.f90, a module containing a few things called by neat.f90, commandline.f90 and output.f90
!(C) Roger Wesson

module mod_functions
use mod_types
use mod_resultarrays
use mod_quicksort

contains

subroutine randomizer(linelist, listlength, norp)
        ! from http://www.netlib.org/random/random.f90
        implicit none
        type(line), dimension(listlength) :: linelist
        integer :: IO, I, j, listlength
        real(kind=dp) :: temp4
        logical, intent(in) :: norp

        real(kind=dp) :: fn_val
        real(kind=dp)     :: s = 0.449871, t = -0.386595, a = 0.19600, b = 0.25472,           &
                    r1 = 0.27597, r2 = 0.27846, u, v, x, y, q
        real(kind=dp) :: half
        real(kind=dp) :: newmean, newsnr, snr

!debugging
#ifdef CO
        print *,"subroutine: randomizer"
#endif

        half = 0.5

        I = 1
        IO=0

        do j = 1,listlength
          do
            call random_number(u)
            call random_number(v)
            v = 1.7156 * (v - half)
            x = u - s
            y = abs(v) - t
            q = x**2 + y*(a*y - b*x)
            if (q .lt. r1) exit
            if (q .gt. r2) cycle
            if (v**2 .lt. -4.0*log(u)*u**2) exit
          enddo
          fn_val = v/u

!                if (j==1) R=3.1+(0.15*fn_val)

          if (linelist(j)%int_err .eq. 0.d0 .and. .not. norp) then ! can't calculate signal to noise for RP calculations
            linelist(j)%intensity = 0.d0
          elseif (linelist(j)%intensity/linelist(j)%int_err .gt. 6.0 .or. norp) then !normal distribution
            temp4=linelist(j)%intensity+(fn_val*linelist(j)%int_err)
            if(temp4 .lt. 0) temp4 = 0.D0
            linelist(j)%intensity = temp4

          elseif (linelist(j)%int_err .ge. linelist(j)%intensity) then !it's an upper limit, take number from semi-gaussian distribution with peak at zero and 5 sigma = intensity
            linelist(j)%intensity = abs(fn_val)*0.2*linelist(j)%intensity
          else !if SN<6, then take lognormal distribution, parameters from Rola & Pelat (1994)
                     !for SN<6, the actual mean is derived from the observed mean using
            snr = linelist(j)%intensity/linelist(j)%int_err
            newmean = 0.0765957/(snr**2) + 1.86037/snr - 0.309695
                     !the actual standard deviation is derived from the observed using
            newsnr = -1.11329/(snr**3) + 1.8542/(snr**2) - 0.288222/snr + 0.18018
                     !(fits to the data in Rola & Pelat's table 6)
                     !the distributions in table 6 give the mean and sigma of log-normal distributions of S/N(obs), given S/N(true).  We don't know S/N(true) but using the distributions as that of the factor by which line fluxes are overestimated is equivalent.  So,
            temp4 = exp(fn_val*newsnr + newmean)
            if (temp4 .lt. 0 ) temp4 = 0.D0
            linelist(j)%intensity = linelist(j)%intensity / temp4

          endif
        enddo

end subroutine randomizer

subroutine init_random_seed()
          integer :: i, n, clock
          integer, dimension(:), allocatable :: seed

!debugging
#ifdef CO
        print *,"subroutine: init_random_seed"
#endif

          n=20
          i=n
          call random_seed(size = n)
          allocate(seed(n))

          call system_clock(count=clock)

          seed = clock + 37 * (/ (i - 1, i = 1, n) /)
          call random_seed(put = seed)

          deallocate(seed)

end subroutine init_random_seed

subroutine write_uncertainties(input_array, uncertainty_array, plaintext, latextext, itemformat, filename, suffix, verbosity,nbins)

!wrapper for the get_uncertainties routine, if called it will do the
!uncertainty calculation and also write the binned and unbinned results to files
!as necessary.
implicit none
real(kind=dp) :: input_array(:)
real(kind=dp), intent(out) :: uncertainty_array(3)
type(arraycount), dimension (:), allocatable :: binned_quantity_result
character(len=35), intent(in) :: plaintext, latextext
character(len=35), intent(in) :: itemformat
character(len=512), intent(in) :: filename
character(len=25), intent(in) :: suffix
logical :: unusual
integer :: verbosity !1 = write summaries, binned and full results; 2=write summaries and binned results; 3=write summaries only
integer :: nbins, i

!debugging
#ifdef CO
        print *,"subroutine: write_uncertainties, quantity=",suffix
#endif

if(size(input_array) .eq. 1) then
!just one iteration, write out the result without uncertainties
  if (input_array(1) .gt. 0.) then
    write (650,itemformat) plaintext,input_array(1)
    write (651,*) latextext," & $", trim(latex_number(input_array(1))),"$\\"
  else
    write (650,*) plaintext,"--"
    write (651,*) latextext," & -- \\"
  endif
else

call get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins)

!write out the binned results with the mode+-uncertainties at the top

if (verbosity .lt. 3) then

  if (allocated(binned_quantity_result) .and. maxval(uncertainty_array) .ne. minval(uncertainty_array)) then
    open(850, file=trim(filename)//"_"//trim(suffix)//"_binned", status='replace',access='sequential', action='write')
    write(unit = 850,fmt=*) uncertainty_array(2),uncertainty_array(2)-uncertainty_array(1),uncertainty_array(3)+uncertainty_array(2)
    write(unit = 850,fmt=*)
  !  do i=1,ii-1
    do i=1,size(binned_quantity_result, 1)
  ! this hacky condition is because for reasons I can't work out right now, the
  ! number of bins allocated to the array is always too large and the last few end
  ! up being full of zeros.  to be fixed soon hopefully.  RW 16/11/2012
!fixed 12/06/2017 - I guess 4.5 years wasn't quite within the hoped for "soon"
!      if (binned_quantity_result(i)%value .gt. 0. .and. binned_quantity_result(i)%counts .gt. 0) then
        write(unit = 850,fmt=*) binned_quantity_result(i)%value,binned_quantity_result(i)%counts
!      endif
    enddo
    close(850)
  endif

endif

!write the unbinned results to file

if (verbosity .lt. 2) then

  open(850, file=trim(filename)//"_"//trim(suffix), status='replace', access='sequential', action='write')
  do i=1,size(input_array)
    write(unit = 850,fmt=*) input_array(i)
  enddo
  close(850)

endif

!write derived value and uncertainties to the summary files

  if (maxval(uncertainty_array) .gt. 0.) then !if this condition is not true, array will be full of zeroes
!    if (unusual) write (650,itemformat) "Warning! Unusual probability distribution.  You should inspect this one:"
    write (650,itemformat) plaintext,uncertainty_array(2),uncertainty_array(3),-uncertainty_array(1)
!    write (651,*) latextext," & $",trim(latex_number(uncertainty_array(2))),"^{",trim(latex_number(uncertainty_array(3))),"}_{",trim(latex_number(-uncertainty_array(1))),"}$ \\"
    if (uncertainty_array(1) .eq. uncertainty_array(3)) then
      write (651,"(A,' & ${',A,'}\pm{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(1)))
    else
      write (651,"(A,' & ${',A,'}^{+',A,'}_{',A,'}$ \\')") latextext,trim(latex_number(uncertainty_array(2))),trim(latex_number(uncertainty_array(3))),trim(latex_number(-uncertainty_array(1)))
    endif
  else
    write (650,*) plaintext,"--  --  --"
    write (651,*) latextext,"& -- & -- & -- \\"
  endif

endif !end of condition checking whether one or more iterations were done

end subroutine write_uncertainties

subroutine get_uncertainties(input_array, binned_quantity_result, uncertainty_array, unusual,nbins)

implicit none
real(kind=dp) :: input_array(:)
real(kind=dp),dimension(:),allocatable :: logarray,exparray
real(kind=dp), intent(out) :: uncertainty_array(3)
real(kind=dp), dimension(:), allocatable :: bintemp
real(kind=dp) :: binsize=0d0
type(arraycount), dimension (:), allocatable, intent(out) :: binned_quantity_result
integer :: arraysize
integer :: nbins
integer :: bincount, i, ii
real(kind=dp) :: comp
real(kind=dp) :: mean=0d0, sd=0d0
real(kind=dp) :: mean_log=0d0, sd_log=0d0
real(kind=dp) :: mean_exp=0d0, sd_exp=0d0
real(kind=dp), dimension(3,3) :: sds=0d0
real(kind=dp) :: tolerance=0d0
real(kind=dp) :: rounding ! value to round to nearest to
logical :: unusual !if not true, then the distribution is normal, log normal or exp normal

!debugging
#ifdef CO
        print *,"subroutine: get_uncertainties. "
#endif

unusual = .false.

if (allocated(binned_quantity_result)) deallocate(binned_quantity_result)

tolerance=0.02 ! to determine whether a probability distribution is normal, log normal or exp normal, we calculate the mean and standard deviation of the original array, the logarithm of the values, and the exponent of the values.  We then determine the fractions of the distribution lying within 1, 2 and 3 standard deviations of the mean in each case.  If these correspond to the 0.68-0.95-0.99 expected for a normal distribution then we know that the distribution is normal, log normal or exp normal and record the uncertainties accordingly.  The fractions must be (0.6827+-tol), (0.9545-tol/2) and (0.9973+-tol/5) for the subroutine to assign a distribution.
sds=0.D0

!binned_quantity_result = 0.D0
uncertainty_array = (/0.0,0.0,0.0/)
arraysize = size(input_array)

!todo: include a verbosity switch
!verbosity: 0=summary file only, 1=summary file + binned results, 2=summary, binned and full
!set verbosity to 0 automatically for single iteration, then we can use the same routine for any n

!sort the input array into ascending order

call qsort(input_array)

! bin the array. use 1/25 x 1 sigma above median as bin size

arraysize = size(input_array)
binsize=input_array(nint(0.841*arraysize))/25d0
nbins=0

if (binsize .gt. 0d0) then

!quantize the input array by taking the integer of each value divided by the binsize, and multiplying by the binsize.

  allocate(bintemp(arraysize))
  bintemp = binsize*nint(input_array/binsize)

  ! count the unique values
  comp=bintemp(1)
  ii=0
  do i=1,arraysize
    if (comp.ne.bintemp(i)) then
      comp=bintemp(i)
      ii=ii+1
    endif
  enddo

  nbins=ii

endif

if (nbins.gt.0) then
  allocate(binned_quantity_result(nbins))
  binned_quantity_result%value=0.D0
  binned_quantity_result%counts=0

!then, go through the quantized array and count the number of occurrences of each value

  comp=bintemp(1)
  bincount=0
  ii=1

  do i=1,arraysize
    if (bintemp(i).ne.comp) then
      binned_quantity_result(ii)%value=comp
      binned_quantity_result(ii)%counts=count(bintemp.eq.comp)
      comp=bintemp(i)
      ii=ii+1
    endif
  enddo

!get median and 1 sigma limits

  uncertainty_array = (/ input_array(nint(0.5d0*arraysize))-input_array(nint(0.159*arraysize)), input_array(nint(0.5d0*arraysize)), input_array(nint(0.841*arraysize))-input_array(nint(0.5d0*arraysize)) /)

!simple test for normality - calculate the mean and standard deviation of the array, determine the number of values within 1, 2 and 3 sigma of the mean.
!if the fractions are close to 68.27, 95.45 and 99.73 then it's normal
!also calculate for log and exp. Temperatures often have an exp-normal distribution but 10**(a typical temperature) exceeds HUGE(0.d0).  Scale to avoid this.

  allocate(logarray(size(input_array)))
  allocate(exparray(size(input_array)))

  where (input_array.gt.0)
    logarray=log10(input_array)
  elsewhere
    logarray=0.d0
  endwhere

  exparray=input_array/maxval(input_array)
  exparray=10**exparray

  mean=sum(input_array)/arraysize
  mean_log=sum(logarray)/arraysize
  mean_exp=sum(exparray)/arraysize

  sd=(sum((input_array-mean)**2)/arraysize)**0.5
  sd_log=(sum((logarray-mean_log)**2,input_array.gt.0.d0)/arraysize)**0.5
  sd_exp=(sum((exparray-mean_exp)**2)/arraysize)**0.5

  sds(1,1)=count(abs(input_array-mean).lt.sd)
  sds(1,2)=count(abs(input_array-mean).lt.2*sd)
  sds(1,3)=count(abs(input_array-mean).lt.3*sd)

  sds(2,1)=count(abs(logarray-mean_log).lt.sd_log)
  sds(2,2)=count(abs(logarray-mean_log).lt.2*sd_log)
  sds(2,3)=count(abs(logarray-mean_log).lt.3*sd_log)

  sds(1,1)=count(abs(exparray-mean_exp).lt.sd_exp)
  sds(1,2)=count(abs(exparray-mean_exp).lt.2*sd_exp)
  sds(1,3)=count(abs(exparray-mean_exp).lt.3*sd_exp)

  sds=sds/arraysize

!scale exp values back up

  mean_exp=mean_exp*maxval(input_array)
  sd_exp=sd_exp*maxval(input_array)

  if (abs(sds(1,1)-0.6827).lt.tolerance .and. abs(sds(1,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(1,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/sd,mean,sd/)
  elseif (abs(sds(2,1)-0.6827).lt.tolerance .and. abs(sds(2,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(2,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/10**(mean_log)-10**(mean_log-sd_log),10**(mean_log),10**(mean_log+sd_log)-10**(mean_log)/)
  elseif (abs(sds(3,1)-0.6827).lt.tolerance .and. abs(sds(3,2)-0.9545).lt.(tolerance*0.5) .and. abs(sds(3,3)-0.9973).lt. (tolerance*0.1)) then
    uncertainty_array(:)=(/log10(mean_exp+sd_exp)-log10(mean_exp),log10(mean_exp),log10(mean_exp-sd_exp)-log10(mean_exp)/)
  else
    unusual = .true.
  endif

  deallocate(bintemp)

else !all results are identical
     !or, almost all of the results are an upper or lower limit with a few
     !actual measurements
     !to account for this second possibility, take the value from the middle of
     !the array as the output quantity.  This will be a useful value in all
     !cases.
  uncertainty_array(1) = 0.D0
  uncertainty_array(2) = input_array(size(input_array)/2)
  uncertainty_array(3) = 0.D0
endif

!round values to 3 significant figures
!if uncertainties get rounded to zero, increment them

if (uncertainty_array(2) .ne. 0.d0 .and. uncertainty_array(2) .ne. 1.d0) then
  if (uncertainty_array(1).ge.0.001 .and.uncertainty_array(1).lt.1.0) then
    rounding=0.001d0
  else
    rounding=10**dble(floor(log10(uncertainty_array(2))))/100.d0
  endif

  uncertainty_array = rounding * nint(uncertainty_array/rounding)
  uncertainty_array(1) = max(rounding,uncertainty_array(1))
  uncertainty_array(3) = max(rounding,uncertainty_array(3))

endif

end subroutine get_uncertainties

character(len=11) function gettime()
implicit none
character(len=10) :: time
character(len=11), save :: oldtime

!debugging
#ifdef CO
        !print *,"function: gettime"
#endif

  call date_and_time(TIME=time)
  gettime = time(1:2)//":"//time(3:4)//":"//time(5:6)//" : "
  if (gettime .eq. oldtime) then
    gettime = "          "
  else
    oldtime = gettime
  endif
  return

end function gettime

character(len=100) function latex_number(inputnumber)
! for any number in the form aEx, write it in latex format, ie a$\times$10$^{x}$
real(kind=dp) :: inputnumber, mantissa
integer :: exponent, pos

!debugging
#ifdef CO
        print *,"function: latex_number, ",inputnumber
#endif

write (latex_number,"(ES14.6)") inputnumber
pos = index (latex_number,'E')
read (latex_number(1:pos-1), '(F6.3)') mantissa
read (latex_number(pos+1:14), '(I3)') exponent

if (exponent .ge. -2 .and. exponent .le. 1) then
  write (latex_number,"(F6.2)") inputnumber !just print out normal number if it's between 0.01 & 10
elseif (exponent .ge. 2 .and. exponent .le. 4) then
  write (latex_number,"(I6)") 10*nint(inputnumber/10) ! write out integer rounded to nearest 10 if it's between 10 and 10,000
else !otherwise, write out a formatted exponent
  write (latex_number,"(F6.2,'\times 10^{',I3,'}')") mantissa,exponent
endif
return

end function latex_number

end module mod_functions
