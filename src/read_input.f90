subroutine readcurrentcycle

	use globalvar
	implicit none
	
	open(1,file ='currentcycle.txt', form = 'formatted', status = 'old')
	read(1,*) icstart
	close(1)

end subroutine readcurrentcycle
subroutine readglobal
! This subroutine is read information from FE_Global.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
		!write(*,*) 'Checking FE_Stations.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Global.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1001, file = 'FE_Global.txt', form = 'formatted', status = 'old')
		read(1001,*) dx
		read(1001,*) sol_op
		read(1001,*) AZTEC_OPTIONS
		read(1001,*) azmaxiter
		read(1001,*) aztol		
		!write(*,*) 'n4nds,nonfs',n4nds, (nonfs(i), i = 1, ntotft), me
	close(1001)
end subroutine readglobal
! #6 readstations --------------------------------------------------------
subroutine readstations1
! This subroutine is read information from FE_Stations.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Stations.txt", EXIST=file_exists)
		!write(*,*) 'Checking FE_Stations.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Stations.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1006, file = 'FE_Stations.txt', form = 'formatted', status = 'old')
		read(1006,*) n4nds
		read(1006,*) (nonfs(i), i = 1, ntotft)
		!write(*,*) 'n4nds,nonfs',n4nds, (nonfs(i), i = 1, ntotft), me
	close(1006)
end subroutine readstations1
! #7 readstations2 --------------------------------------------------------
subroutine readstations2
! This subroutine is read information from FE_Stations.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="FE_Stations.txt", EXIST=file_exists)
		!write(*,*) 'Checking FE_Stations.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="FE_Stations.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'FE_Stations.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1006, file = 'FE_Stations.txt', form = 'formatted', status = 'old')
		read(1006,*) 
		read(1006,*) 
		read(1006,*)
		do i = 1, ntotft
			do j = 1, nonfs(i)
				read(1006,*) xonfs(1,j,i), xonfs(2,j,i)
				!write(*,*) xonfs(1,j,i), xonfs(2,j,i)
			enddo 
			read(1006,*)
		enddo
		do i = 1, n4nds
			read(1006,*) x4nds(1,i), x4nds(2,i), x4nds(3,i)
			!write(*,*) x4nds(1,i), x4nds(2,i), x4nds(3,i)
		enddo 
	close(1006)
	
	xonfs=xonfs*1000.0d0  !convert from km to m
	x4nds=x4nds*1000.0d0		
		
end subroutine readstations2