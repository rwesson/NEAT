program makedata

implicit none

double precision, dimension(46) :: temp
double precision, dimension(22) :: wavelengths
character*7, dimension(5) :: filenames
integer :: io, i, j, k
character*2 :: region
character*9 :: outfile
integer, dimension(5) :: lengths

wavelengths = (/ 3727.00, 3868.75, 4101.74, 4072.00, 4340.47, 4471.50, 4958.91, 5006.84, 5200.00, 5754.60, 5875.66, 6300.30, 6312.10, 6548.10, 6562.77, 6583.50, 6678.16, 6716.44, 6730.82, 7135.80, 7325.00, 9068.60 /)
filenames = (/"ngc1232","ngc1365","ngc2903","ngc2997","ngc5236"/)
lengths = (/ 13,15,9,14,18 /)

do i=1,5
  open (100, file=trim(filenames(i)), iostat=io)
  do j=1,lengths(i)
    read (100,*) temp
    write (region,"(I2.2)") nint(temp(2))
    open (101, file=trim(filenames(i))//"_"//region//".dat")
    do k=1,22
      write (101,"(3(F7.2))") wavelengths(k),temp(k*2+1),temp(k*2+2)
    end do
    write (101,"(3(F7.2))") 4861.33,100.0,2.0
  end do
end do

end program makedata
