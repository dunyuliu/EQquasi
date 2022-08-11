subroutine readcurrentcycle

	use globalvar
	implicit none
	
	open(1,file ='currentcycle.txt', form = 'formatted', status = 'old')
	read(1,*) icstart
	close(1)

end subroutine readcurrentcycle

! #2 readmodelgeometry -------------------------------------------------
subroutine readmodel
! This subroutine is read information from FE_global.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	
	if (me == 0) then 
		INQUIRE(FILE="model.txt", EXIST=file_exists)
		!write(*,*) 'Checking FE_Model_Geometry.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'model.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="model.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'model.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1002, file = 'model.txt', form = 'formatted', status = 'old')
		read(1002,*) xmin, xmax
		read(1002,*) ymin, ymax
		read(1002,*) zmin, zmax
		read(1002,*) xminc, xmaxc, zminc ! creeping zone boundaries.
		read(1002,*) dis4uniF, dis4uniB
		read(1002,*) rat
		read(1002,*) dx 
		read(1002,*) mat0(1,1), mat0(1,2), mat0(1,3)
		read(1002,*) C_elastic
		read(1002,*) friclaw
		read(1002,*) ntotft ! 10
		read(1002,*) bp
		read(1002,*) ksi, minDc
		read(1002,*) slipr_thres
		read(1002,*) far_load_rate
		read(1002,*) load_slip_rate
		read(1002,*) init_norm
		read(1002,*) sol_op
		read(1002,*) AZTEC_OPTIONS
		read(1002,*) azmaxiter
		read(1002,*) aztol
	close(1002)
	
	if (xminc < xmin .or. xmaxc > xmax .or. zminc < zmin) stop ! Creeping zone bounds should be within model bounds. 
	!ccosphi=coheplas*dcos(atan(bulk))
	!sinphi=dsin(atan(bulk))
	!nstep=idnint(term/dt)
	!rdampk=rdampk*dt	
	!tv = 2.0d0*dx/3464.0d0
	
end subroutine readmodel

! #5 readfric --------------------------------------------------------
subroutine readfric
! This subroutine is read information from FE_Material.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="fric.txt", EXIST=file_exists)
		!write(*,*) 'Checking FE_Fric.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'fric.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="fric.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'fric.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1005, file = 'fric.txt', form = 'formatted', status = 'old')
	if (friclaw == 1) then ! slip-weakening
		read(1005,*) fric_sw_fs
		read(1005,*) fric_sw_fd
		read(1005,*) fric_sw_D0
	endif 
	if (friclaw == 3) then ! aging law
		do i = 1, 4
			read(1005,*)
		enddo
		read(1005,*) fric_rsf_a
		read(1005,*) fric_rsf_deltaa0
		read(1005,*) fric_rsf_b
		read(1005,*) fric_rsf_Dc
		read(1005,*) fric_rsf_r0
		read(1005,*) fric_rsf_v0
		read(1005,*)
		read(1005,*) fric_rsf_vinix	
		read(1005,*) fric_rsf_viniz	
	endif 	
	if (friclaw == 4) then ! strong-rate weakening
		do i = 1, 4
			read(1005,*)
		enddo
		read(1005,*) fric_rsf_a
		read(1005,*) fric_rsf_deltaa0
		read(1005,*) fric_rsf_b
		read(1005,*) fric_rsf_Dc
		read(1005,*) fric_rsf_r0
		read(1005,*) fric_rsf_v0
		read(1005,*)
		read(1005,*) fric_rsf_vinix	
		read(1005,*) fric_rsf_viniz	
		read(1005,*)
		read(1005,*) fric_rsf_fw
		read(1005,*) fric_rsf_vw
		read(1005,*) fric_rsf_deltavw0	
	endif 		
	if (friclaw == 5) then ! strong-rate weakening + termop
		do i = 1, 4
			read(1005,*)
		enddo
		read(1005,*) fric_rsf_a
		read(1005,*) fric_rsf_deltaa0
		read(1005,*) fric_rsf_b
		read(1005,*) fric_rsf_Dc
		read(1005,*) fric_rsf_r0
		read(1005,*) fric_rsf_v0
		read(1005,*)
		read(1005,*) fric_rsf_vinix	
		read(1005,*) fric_rsf_viniz	
		read(1005,*)
		read(1005,*) fric_rsf_fw
		read(1005,*) fric_rsf_vw
		read(1005,*) fric_rsf_deltavw0		
		read(1005,*)
		read(1005,*) fric_tp_a_th
		read(1005,*) fric_tp_rouc
		read(1005,*) fric_tp_lambda
		read(1005,*) fric_tp_h
		read(1005,*) fric_tp_a_hy
		read(1005,*) fric_tp_deltaa_hy0
		read(1005,*) fric_ww
		read(1005,*) fric_w
		read(1005,*) 
		read(1005,*) fric_ini_sliprate
		read(1005,*) fric_tp_Tini
		read(1005,*) fric_tp_pini
		! read(1005,*) dxtp
		! read(1005,*) tpw
	endif 	
	close(1005)
end subroutine readfric
! #6 readstations --------------------------------------------------------
subroutine readstations1
! This subroutine is read information from stations.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="stations.txt", EXIST=file_exists)
		!write(*,*) 'Checking stations.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'stations.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="stations.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'stations.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1006, file = 'stations.txt', form = 'formatted', status = 'old')
		read(1006,*) n4nds
		read(1006,*) (nonfs(i), i = 1, ntotft)
		!write(*,*) 'n4nds,nonfs',n4nds, (nonfs(i), i = 1, ntotft), me
	close(1006)
end subroutine readstations1
! #7 readstations2 --------------------------------------------------------
subroutine readstations2
! This subroutine is read information from stations.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j 
	
	if (me == 0) then 
		INQUIRE(FILE="stations.txt", EXIST=file_exists)
		!write(*,*) 'Checking stations.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'stations.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="stations.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'stations.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1006, file = 'stations.txt', form = 'formatted', status = 'old')
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

! #8 read_rough_geometry ------------------------------------------------
subroutine read_fault_rough_geometry
! This subroutine is read information from FE_Fault_Rough_Geometry.txt
	use globalvar
	implicit none
	include 'mpif.h'

	logical::file_exists
	integer(kind=4):: i, j
	
	if (me == 0) then 
		INQUIRE(FILE="rough_geo.txt", EXIST=file_exists)
		!write(*,*) 'Checking rough_geo.txt by the master procs', me
		if (file_exists == 0) then
			write(*,*) 'rough_geo.txt is required but missing ...'
			pause
		endif 
	endif 
	if (me == 0) then 
		INQUIRE(FILE="rough_geo.txt", EXIST=file_exists)
		if (file_exists == 0) then
			write(*,*) 'rough_geo.txt is still missing, so exiting EQdyna'
			stop
		endif 
	endif 	
	
	open(unit = 1008, file = 'rough_geo.txt', form = 'formatted', status = 'old')
		read(1008,*) nnx, nnz
		read(1008,*) dxtmp 
	close(1008)
	
	allocate(rough_geo(3,nnx*nnz))
	
	open(unit = 1008, file = 'rough_geo.txt', form = 'formatted', status = 'old')
		read(1008,*)
		read(1008,*)
		do i = 1, nnx*nnz
			read(1008,*) rough_geo(1,i), rough_geo(2,i), rough_geo(3,i)
		enddo
	close(1008)	
		
end subroutine read_fault_rough_geometry