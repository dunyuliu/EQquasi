!===================================================================!
program eqquasi3d
use globalvar
implicit none
include 'mpif.h'
real(kind=8)::timebegin,timeover
integer(kind=4)::numnp,numel,neq,IERR,me,nprocs
real(kind=8),allocatable,dimension(:,:)::x,mat
real(kind=8),allocatable,dimension(:)::brhs,d,v
integer(kind=4),allocatable,dimension(:,:)::ien,id
integer(kind=4),allocatable,dimension(:)::et 
!For faulting
integer(kind=4),allocatable,dimension(:)::nftnd
integer(kind=4),allocatable,dimension(:,:) :: anonfs
integer(kind=4),allocatable,dimension(:,:,:) :: nsmp
real(kind=8),allocatable,dimension(:,:) :: fnft,arn,r4nuc,arn4m,slp4fri
real(kind=8),allocatable,dimension(:,:,:) :: fric,un,us,ud,fltslp
integer(kind=4)::nftmx,nonmx,n4onf
integer(kind=4),dimension(1)::nonfs=(/9/)
real(kind=8),dimension(2,9,1)::xonfs

!integer istatus(MPI_STATUS_SIZE)

CALL MPI_INIT(IERR)
call mpi_comm_rank(MPI_COMM_WORLD,me,IERR)
call mpi_comm_size(MPI_COMM_WORLD,nprocs,IERR)

fltxyz(1,1,1)=-30.0d3
fltxyz(2,1,1)=30.0d3
fltxyz(1,2,1)=0.0d0
fltxyz(2,2,1)=0.0d0
fltxyz(1,3,1)=-30.0d3
fltxyz(2,3,1)=0.0d0
fltxyz(1,4,1)=270.0d0/180.0d0*pi
fltxyz(2,4,1)=90.0d0/180.0d0*pi
!nstep=int(term/dt)	
!...time hostories output steps
nplpts = 0	!initialize number of time history plot
if (nhplt > 0) then
	nplpts = int(nstep/nhplt) + 2
endif
!======================1=========================!
!================MESHGENRATION===================!
allocate(nftnd(ntotft))
if (me==0) write(*,*) 'BEGORE MESH4NUM'
	call mesh4num(numnp,numel,neq,nftnd,me)
if (me==0) write(*,*) 'AFTER MESH4NUM AND INITIATE PARAMETERS FOR FEM'

allocate(x(ndof,numnp),id(ndof,numnp),ien(8,numel),mat(numel,5),et(numel))

nftmx=maxval(nftnd) !max fault nodel num for all faults, used for arrays.
if(nftmx<=0) nftmx=1  !fortran arrays cannot be zero size,use 1 for 0
nonmx=sum(nonfs)    !max possible on-fault stations number
allocate(nsmp(2,nftmx,ntotft),fnft(nftmx,ntotft),un(3,nftmx,ntotft),&
			us(3,nftmx,ntotft),ud(3,nftmx,ntotft),fric(50,nftmx,ntotft),&
			arn(nftmx,ntotft),r4nuc(nftmx,ntotft),anonfs(3,nonmx),&
			slp4fri(nftmx,ntotft),fltslp(3,nftmx,ntotft))
fnft=10000.d0!Should be initialized over 600.
fric = 0.0d0
un = 0.0d0
us = 1000.0d0
ud = 0.0d0
arn = 0.0d0
r4nuc = 0.0d0
anonfs = 0
slp4fri = 0.0d0
fltslp = 0.0d0

if (me==0) write(*,*) 'AFTER INITIATION AND BEFORE MESHGEN'	
call meshgen(numnp,numel,neq,ien,mat,et,x,id,&
			nftnd,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,me)
if (me==0) write(*,*) 'AFTER MESHGEN AND BEFORE CRSELEM'	

!======================2=========================!
!================STIFFNESS COMP==================!	
call CRSelem(numel,numnp,neq,x,mat,ien,id,&
			nftnd,n4onf,xonfs,nonfs,&
			nftmx,nonmx,nsmp,un,us,ud,fric,arn,r4nuc,&
			anonfs,fnft,slp4fri,et,me)
if (me==0) write(*,*) 'AFTER CRSELEM AND THE END OF THE SIMULATION'

!======================3=========================!
!=================OUTPUT PHASE===================!
!CALL MPI_FINALIZE(IERR)
STOP
end 