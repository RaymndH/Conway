program row
use mpi
integer numtasks, taskid, len, ierr, status(MPI_STATUS_SIZE)
character(MPI_MAX_PROCESSOR_NAME) hostname
character(len=100) :: filename
integer msg
integer m,n,s,i,j,k,si,REALSIZE
integer, dimension(2) :: sizes,subsizes,starts
integer, allocatable :: widths(:), scounts(:), displs(:),prints(:) 
real, allocatable :: A(:,:), B(:,:), buffer(:,:), transbuffer(:,:),lilbuffer(:)
!write(*,*) 'mpicommworld is', MPI_COMM_WORLD
integer :: otype, rowtype, recvtype, newtype, resizedtype
integer(kind=MPI_ADDRESS_KIND) :: start, extent
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
call MPI_GET_PROCESSOR_NAME(hostname,len,ierr)
 ! write(*,*)'program started'
m=20
n=20
allocate(prints(3))
prints=(/20,40,80/)
if (taskid==0) then
  write(*,*) 'Printing at ',prints
   write(*,*) 'Initial State:'
  allocate(B(m,n))
! Random configuration
  call random_number(B)
  where (B<.5) B=0
  where (B>.5) B=1
! ordered grid
  k=0
  do i=1,m
    do j=1,n
      B(i,j)=k
      k=k+1
    end do
  end do
! beacon
!  B=0
!  B(2:3,2:3)=1
!  B(4:5,4:5)=1
! Glider
  B=0
  B(2,1)=1
  B(3,2)=1
  B(1:3,3)=1
  do i=1,m
    write(*,'(999(f4.0))') (B(i,j), j=1,n)
  end do
!Split up into the vectors
end if
  
  s=ceiling(real(m)/numtasks)
  allocate(widths(numtasks),scounts(numtasks),displs(numtasks))
  widths=s
  widths(numtasks)= m-s*(numtasks-1)
  scounts=widths*n
  displs(1)=0
  do i=2,numtasks
    displs(i)=displs(i-1)+widths(i-1)
  end do

!if (taskid==0) then
if (1==0) then
  write(*,*) ''
  write(*,*) 'widths', widths
  write(*,*) 'displs',  displs
endif
  !msg=taskid
  si=widths(taskid+1)
  allocate(transbuffer(n,si))
  allocate(buffer(si,n))
  buffer=7
  allocate(A(si+2,n+2))
  A=0
  allocate(lilbuffer(n))
  lilbuffer=0
!------ Create derived datatypes
  !call MPI_Type_vector(n,1,si,MPI_real,otype,ierr)
  call MPI_Type_vector(3,1,2,MPI_real,rowtype,ierr)
  call MPI_Type_SIZE(MPI_REAL,REALSIZE,IERR)
  extent = 1*REALSIZE
  !write(*,*) extent
  start=0
  !call MPI_Type_create_resized(otype,start,extent,rowtype,ierr)
  call MPI_Type_commit(rowtype, ierr)
 ! call MPI_Type_vector(n,1,si,MPI_real,otype,ierr)

  !call MPI_Type_create_resized(otype,start,extent,recvtype,ierr)
 
!---------
  sizes=(/m,n/)
  starts=(/0,0/)
  subsizes=(/si,n/)
  call MPI_Type_create_subarray(2,sizes,subsizes,starts,&
        MPI_ORDER_FORTRAN,MPI_REAL,newtype,ierr)
  call MPI_Type_size(MPI_REAL,REALSIZE,ierr)
  extent=si*REALSIZE
  extent=1*REALSIZE
 ! write(*,*)'extent', extent
  call MPI_Type_create_resized(newtype,start,extent,resizedtype,ierr)
  call MPI_Type_commit(resizedtype, ierr)
!----------
  widths=1 
!----------  
  !call MPI_TYPE_VECTOR()
 
 ! write(*,*) 'calling scatterv pls'
!!------ Initial SCATTERV
  call MPI_SCATTERv(B,widths,displs,resizedtype,         &
    buffer,si*n,MPI_REAL,0,MPI_COMM_WORLD,ierr)
  !write(*,*) 'Server', taskid
!------ do loop which runs entire system kmax times
  if (0==1) then
   write(*,*) 'initial buffer on server',taskid
   !buffer=transpose(transbuffer)
   do i=1,si
    write(*,'(999(f4.0))') (buffer(i,j), j=1,n)
   enddo
  write(*,*)'all done' 
  endif
   
do k=1,prints(size(prints))
!------ Pictureframing step
 ! write(*,*) 'About to start pictureframing step ',k
  !write(*,*) 'Size of A is', size(A)
  ! write(*,*) 'A iter ',k,'server',taskid
  !  do i=1,m+2
  !  write(*,'(999(f4.0))') (A(i,j), j=1,si+2)
  !  enddo 

  !write(*,*) 'A w/out side boundaries, iter',k,'on server',taskid
  !do i=1,m+2
  !  write(*,'(999(f4.0))') (A(i,j), j=1,si+2)
  !enddo
  !write(*,*) 'A end'   
  !write(*,*) 'Preparing for sendrecv on',taskid,'at iteration',k
  ! write(*,*) 'buffer w/out side boundaries, iter',k,'on server',taskid
  do i=1,si
    !write(*,'(999(f4.0))') (buffer(i,j), j=1,n)
  enddo
  if (1==0) then
  ! write(*,*) 'please god lol'
  ! write(*,*) buffer(1,:)
  ! write(*,*) 'yeah lol'
  endif
  A(2:si+1,2:n+1)=buffer
  lilbuffer=buffer(1,:) 
  !A(2:si+1,2:n+1)=buffer
  !write(*,*)'lilbuffer on server',taskid,' is '
  !write(*,*)  lilbuffer
  !write(*,*) 'end lilbuffer'
   call mpi_sendrecv(lilbuffer,n,mpi_real,modulo(taskid+1,numtasks),0,&
                    lilbuffer,n,mpi_real,modulo(taskid-1,numtasks),0,MPI_COMM_WORLD,status,ierr)
  A(si+2,2:n+1)=lilbuffer
  !write(*,*)'lilbuffer on server',taskid,' is '
  !write(*,*)  lilbuffer
  !write(*,*) 'end lilbuffer'
 lilbuffer=buffer(si,:)
  call mpi_sendrecv(lilbuffer,n,mpi_real,modulo(taskid-1,numtasks),0,&
                    lilbuffer,n,mpi_real,modulo(taskid+1,numtasks),0,MPI_COMM_WORLD,status,ierr)
  A(1,2:n+1)=lilbuffer  
 A(:,n+2)=A(:,2)
  A(:,1)=A(:,n+1)
   !write(*,*)'A on',taskid,' is '
  do i=1,si+2
    !write(*,'(99(f4.0))')  (A(i,j), j=1,n+2)
  enddo
    if (0==1) then
  write(*,*) 'A with all boundary is:'
  do i=1,si+2
    write(*,'(99(f4.0))')  (A(i,j), j=1,n+2)
  enddo
  endif
 !write(*,*) 'sendrecv on',taskid,'at iteration',k,' worked'
 !write(*,*) 'A with side ghosts, iteration', k
 ! do i=1,m+2
 !   write(*,'(99(f4.0))') (A(i,j), j=1,si+2)
 ! end do
  !write(*,*) ''
! Game of life loop
!!!----- Uncomment this when you want it to actually work!!!
   do i=1,si
     do j=1,n
       buffer(i,j)=sum(A(i:i+2,j:j+2))-A(i+1,j+1)
     end do
   end do
   where((buffer.ne.2).and.(buffer.ne.3) ) A(2:si+1,2:n+1)=0
   where(Buffer.eq.3) A(2:si+1,2:n+1)=1     
   Buffer=A(2:si+1,2:n+1)  
!!!----- End of uncomment
  !write(*,*) 'server ',taskid,' A after GoL step ',k,' is'
  !do i=1,m+2
  !  write(*,'(99(f4.0))') (A(i,j), j=1,si+2)
  !enddo
!------ block which prints to output
!  A=A+1
!   if (1==0) then
  if (any(prints==k)) then
   call MPI_GATHERV(Buffer,si*n,MPI_Real,         &
    B,widths,displs,resizedtype,0,MPI_COMM_WORLD,ierr)
  if(taskid==0) then
    write(*,*) 'iteration',k
    write(*,*) ''
    do i=1,m
      write(*,'(999(f4.0))') (B(i,j), j=1,n)
    enddo
  end if
  endif

  !write(*,*) 'A on taskid',taskid,'k=',k
  !do i=1,m+2
  !  write(*,'(999(f4.0))') (A(i,j),j=1,si+2)
  !enddo
  !write(*,*) ''
  !write(*,*) 'about to finish iteration ',k
!---- end print output
end do
!------ Everything here up must be iterated in the game of life loop for k
!------ iterations
  deallocate(buffer)
if(taskid==0) then
  !write(*,*) s
  deallocate(B)
end if
call MPI_FINALIZE(ierr)
end program

