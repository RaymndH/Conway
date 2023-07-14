program Column_Decomp
use mpi
integer numtasks, taskid, len, ierr, status(MPI_STATUS_SIZE)
character(MPI_MAX_PROCESSOR_NAME) hostname
character(len=100) :: filename
integer msg
integer m,n,s,i,j,k,si
integer, allocatable :: widths(:), scounts(:), displs(:),prints(:) 
real, allocatable :: A(:,:), B(:,:), buffer(:,:), lilbuff(:)
!write(*,*) 'mpicommworld is', MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,taskid,ierr)
call MPI_GET_PROCESSOR_NAME(hostname,len,ierr)

m=20
n=20
allocate(prints(3))

prints=(/20,40,80/)

!allocate(buffer(m,n))
if (taskid==0) then
  write(*,*) 'system prints at', prints
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
  B=0
  B(2:3,2:3)=1
  B(4:5,4:5)=1
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
  s=ceiling(real(n)/numtasks)
  allocate(widths(numtasks),scounts(numtasks),displs(numtasks))
  widths=s
  widths(numtasks)= n-s*(numtasks-1)
  scounts=widths*m
  displs(1)=0
  do i=2,numtasks
    displs(i)=displs(i-1)+scounts(i-1)
  end do
!if (taskid==0) then
if (1==0) then
  write(*,*) ''
  write(*,*) 'scounts', scounts
  write(*,*) 'displs',  displs
endif
  !msg=taskid
  si=widths(taskid+1)
  allocate(buffer(m,si))
  buffer=0
  allocate(A(m+2,si+2))
  A=0
!------ Initial SCATTERV
  call MPI_SCATTERv(B,scounts,displs,MPI_FLOAT,         &
    A(2:m+1,2:si+1),scounts(taskid+1),MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
!  write(*,*) 'Server', taskid
!------ do loop which runs entire system kmax times
do k=1,prints(size(prints))
!------ Pictureframing step
  !write(*,*) 'About to start pictureframing step ',k
  !write(*,*) 'Size of A is', size(A)
  ! write(*,*) 'A iter ',k,'server',taskid
  !  do i=1,m+2
  !  write(*,'(999(f4.0))') (A(i,j), j=1,si+2)
  !  enddo 
  A(m+2,:)=A(2,:)
  A(1,:)=A(m+1,:)
  !write(*,*) 'A w/out side boundaries, iter',k,'on server',taskid
  !do i=1,m+2
  !  write(*,'(999(f4.0))') (A(i,j), j=1,si+2)
  !enddo
  !write(*,*) 'A end'   
!  write(*,*) 'Preparing for sendrecv on',taskid,'at iteration',k
  call mpi_sendrecv(A(:,2),m+2,mpi_real,modulo(taskid-1,numtasks),0,&
                    A(:,si+2),m+2,mpi_real,modulo(taskid+1,numtasks),0,MPI_COMM_WORLD,status,ierr)
  call mpi_sendrecv(A(:,si+1),m+2,mpi_real,modulo(taskid-1,numtasks),0,&
                    A(:,1),m+2,mpi_real,modulo(taskid+1,numtasks),0,MPI_COMM_WORLD,status,ierr)
! write(*,*) 'sendrecv on',taskid,'at iteration',k,' worked'
 !write(*,*) 'A with side ghosts, iteration', k
 ! do i=1,m+2
 !   write(*,'(99(f4.0))') (A(i,j), j=1,si+2)
 ! end do
 ! write(*,*) ''
! Game of life loop
!!!----- Uncomment this when you want it to actually work!!!
  do i=1,m
    do j=1,si
      buffer(i,j)=sum(A(i:i+2,j:j+2))-A(i+1,j+1)
    end do
  end do
  where((buffer.ne.2).and.(buffer.ne.3) ) A(2:m+1,2:si+1)=0
  where(Buffer.eq.3) A(2:m+1,2:si+1)=1     
  Buffer=A(2:m+1,2:si+1)  
!!!----- End of uncomment
  !write(*,*) 'server ',taskid,' A after GoL step ',k,' is'
  !do i=1,m+2
  !  write(*,'(99(f4.0))') (A(i,j), j=1,si+2)
  !enddo
!------ block which prints to output
!  A=A+1
!   if (0==0) then
if (1==0) then
  write(*,*) k
endif
 if(any(prints==k)) then
  !write(*,*)'k=',k
   call MPI_GATHERV(A(2:m+1,2:si+1),scounts(taskid+1),MPI_FLOAT,         &
    B,scounts,displs,MPI_FLOAT,0,MPI_COMM_WORLD,ierr)
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
subroutine printmat(A,m,n)
    implicit none
    integer :: m,n
    real :: A(m,n)
    integer i,j
    write(*,*) 'm is ', m
    write(*,*) 'n is ', n
    write(*,*) A
    write(*,*) 'end of matrix'
 do i=1,m
    write(*,*) A(i,1)
    write(*,'(999(f4.0))') (A(i,j), j=1,n)
  end do
end subroutine
subroutine test()
  write(*,*) 'subroutine test called'
end subroutine test
