program main
  implicit none
  character(len=100) :: info1,info2
  integer*4  :: data(5)
  integer*4 :: rec_tot,npts
  integer*4 :: i,read_num,remainder
  real*8,dimension(:),allocatable  :: area
  real*8,dimension(:),allocatable  :: coef
  integer*4,dimension(:),allocatable :: keep_int_info1,keep_int_info2,keep_int_info3
!  open(unit=15,file='./data/small_scale/unstructured_test.stor')
!  open(unit=15,file='./data/structured_large_scale/structured_10th.stor')
!  open(unit=15,file='./data/structured_large_scale/structured_1mln.stor')
!  open(unit=15,file='./data/unstructured_large_scale/unstructured_10th.stor')
!  open(unit=15,file='./data/unstructured_large_scale/unstructured_1mln.stor')
!  open(unit=15,file='./data/4_fracture_network/fract4_dfn.stor')
!  open(unit=15,file='./data/2_structured_fracture_network/structured_intersection.stor')
!  open(unit=15,file='full.stor')


!  open(unit=15,file='../../data/small_scale/structured_triplane.stor')
!  open(unit=15,file='../../../data/DFN_networks/2_structured_fracture_network/structured_intersection.stor')
!  open(unit=15,file='../../../data/DFN_networks/3_tri_fracture_network/full_mesh.stor')
!  open(unit=15,file='../../../data/DFN_networks/4_fractures/full_mesh.stor')
  open(unit=15,file='../../data/DFN_networks/100_fractures/full_mesh.stor')
!   open(unit=15,file='../../data/DFN_networks/1000_fractures/full_mesh.stor')
!   open(unit=15,file='../../data/DFN_networks/17237_fractures/full_mesh.stor')


!17 format(5(1PE20.12))
17 format(5(F12.4))
18 format(5i15)
  read (15,'(a)') info1
  read (15,'(a)') info2
  read (15,*) data
  rec_tot = int(data(1))
  npts  = int(data(2))
  open(unit=37,file='test.stor')            
  print*, rec_tot,npts  

  allocate(area(npts))
  print*,"done 41"
  allocate(keep_int_info1(npts+1+rec_tot))
  allocate(keep_int_info2(npts+1+rec_tot))
  print*, "done 43"
  allocate(keep_int_info3(npts))
  print*, "done 45"
  allocate(coef(rec_tot))
  print*,"done 47"
  
  read_num = int(npts/5)
  remainder = npts-read_num*5
  do i=1,read_num
    read(15,"(5F20.12)") area((i-1)*5+1:(i-1)*5+5)
  enddo
  if (remainder>0)   read(15,*) area(read_num*5+1:npts)
  
  print*,"area",area(npts-9:npts)

  read_num = int((npts+1+rec_tot)/5)
  remainder = npts+1+rec_tot - read_num*5
  do i=1,read_num
    read(15,"(5I10)") keep_int_info1((i-1)*5+1:(i-1)*5+5)
  enddo
  if (remainder>0)  read(15,*) keep_int_info1(read_num*5+1:npts+1+rec_tot)
  print*,"info1",keep_int_info1(npts+1+rec_tot-9:npts+1+rec_tot)

  read_num = int((npts+1+rec_tot)/5)
  remainder = npts+1+rec_tot-read_num*5
  do i=1,read_num
    read(15,"(5I10)") keep_int_info2((i-1)*5+1:(i-1)*5+5)
  enddo
  if (remainder>0)  read(15,*) keep_int_info2(read_num*5+1:npts+1+rec_tot)
  print*,"info2",keep_int_info2(npts+1+rec_tot-9:npts+1+rec_tot)

  read_num = int(npts/5)
  remainder = npts-read_num*5
  do i=1,read_num
    read(15,"(5I10)") keep_int_info3((i-1)*5+1:(i-1)*5+5)
  enddo
  if (remainder>0) read(15,*) keep_int_info3(read_num*5+1:npts)
  print*,"info3",keep_int_info3(npts-9:npts)
   
  read_num = int(rec_tot/5)
  remainder = rec_tot-read_num*5
  do i=1,read_num
    read(15,"(5F20.12)") coef((i-1)*5+1:(i-1)*5+5)
  enddo
  if (remainder>0)  read(15,*) coef(read_num*5+1:rec_tot)
  print*,"coef",coef(rec_tot-9:rec_tot)
  close(15)


!  read (15,*) keep_int_info1
!  print*, "done 54"
!  read (15,*) keep_int_info2
!  print*, "done 56"
!  read (15,*) keep_int_info3
!  print*, "done 58"
!  read (15,*) coef
!  print*, "done 60"
  
  write(37,'(a)') info1
  write(37,'(a)') fdate()
  write(37,18) data       
  write(37,17),(dabs(area(i)),i=1,npts)
  write(37,18),(keep_int_info1(i),i=1,npts+1+rec_tot)
  write(37,18),(keep_int_info2(i),i=1,npts+1+rec_tot)
  write(37,18),(keep_int_info3(i),i=1,npts)
  write(37,17) (dabs(coef(i)),i=1,rec_tot)
  close(37)
end program main
