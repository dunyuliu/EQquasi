! netcdf_io reads and writes out arrays into netcdf format.
! A Fortran example code is shown here https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_wr.f90
! A companion exmaple code to read the data written by sfc_pres_temp_wr.f90 could be found via 
! https://www.unidata.ucar.edu/software/netcdf/examples/programs/sfc_pres_temp_rd.f90.
! Note that if a 6 X 12 data on a lat-lon grid is to be created, the variable array should be in 
! this format - on_fault_vars(lon, lat). 
! In our cases, lat = z/dip and lon = x/strike. So, on_fault_vars should be on_fault_vars(nx,nz)

! subroutines contained in this file includes: 
! - #1: netcdf_write
! - #2: netcdf_write_on_fault
! - #3: netcdf_write_roughness
! - #4: netcdf_read_on_fault
! - #5: netcdf_read_on_fault_restart
! - #A1: check

! #1
subroutine netcdf_write(outfile, outtype)
    use netcdf
    use globalvar
    implicit none 
    character (len = 50 ) :: outfile, outtype, lat_name, lon_name, lat_units, lon_units, var_1_name, var_1_unit, UNITS
    integer (kind = 4) :: ncid, lat_dimid, lon_dimid, lat_varid, lon_varid, var_1_id, ilat, ilon, i, j, nlat, nlon
    integer (kind = 4) :: dimids(2),  nsmp_tmp(2,nftmx)
    integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index

    UNITS = 'units'
    if (outtype == 'disp') then 
        lat_name = 'node_id'
        lon_name = 'disp_dim_id'
        lat_units = 'unit'
        lon_units = 'unit'
        var_1_name = 'disp'
        var_1_unit = 'meters'
        nlat = numnp
        nlon = ndof
    elseif (outtype == 'coor') then 
        lat_name = 'node_id'
        lon_name = 'coor_dim_id'
        lat_units = 'unit'
        lon_units = 'unit'
        var_1_name = 'coordinates'
        var_1_unit = 'meters'
        nlat = numnp
        nlon = ndof
    elseif (outtype == 'ien') then 
        lat_name = 'element_id'
        lon_name = 'node_index_id'
        lat_units = 'unit'
        lon_units = 'unit'
        var_1_name = 'node_id'
        var_1_unit = 'unit'
        nlat = numel
        nlon = nen
    elseif (outtype == 'nsmp') then
        lat_name = 'fault_node_id'
        lon_name = 'slave_master_id'
        lat_units = 'unit'
        lon_units = 'unit'
        var_1_name = 'node_id'
        var_1_unit = 'unit'
        nlat = nftmx
        nlon = 2
        do i = 1, nftmx
            do j = 1, 2
                nsmp_tmp(j,i) = nsmp(j,i,1)
            enddo 
        enddo 
    endif 
    allocate(lat_index(nlat))
    allocate(lon_index(nlon))
    lat_index = (/ (i, i = 1, nlat) /)
    lon_index = (/ (i, i = 1, nlon) /)
    ! Create the netCDF file.
    call check(nf90_create(outfile, NF90_CLOBBER, ncid))
    
    ! Define the dimensions.
    call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
    call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
    
    ! Define coordiante variables. They will hold the coordinate 
    ! information, that is, the latitudes (y), and longitudes (x). A varid is 
    ! returned for each.
    call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
    call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
    
    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
    call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )    

    ! Define the netcdf variables. The dimids array is used to pass the 
    ! dimids of the dimensions of the netCDF variables.
    dimids = (/ lon_dimid, lat_dimid /)
    if (outtype == 'disp' .or. outtype == 'coor') then 
        call check(nf90_def_var(ncid, var_1_name, NF90_REAL, dimids, var_1_id)) 
    elseif (outtype == 'nsmp' .or. outtype == 'ien') then 
        call check(nf90_def_var(ncid, var_1_name, NF90_INT, dimids, var_1_id)) 
    endif 
    ! Assign units attributes to the pressure and temperature netCDF
    ! variables.
    call check(nf90_put_att(ncid, var_1_id, UNITS, var_1_unit))
    
    ! End definitions.
    call check(nf90_enddef(ncid))
    ! Write data.
    ! Write the coordinate variable data. This will put the x, and y 
    ! of our data grid into the netCDF file.
    call check(nf90_put_var(ncid, lat_varid, lat_index))
    call check(nf90_put_var(ncid, lon_varid, lon_index))
    
    ! Write the data. This will write our displacement fields "cons" which is defined 
    ! globally. Its dimension is 3 by numnp (row by column).
    if (outtype == 'disp') then 
        call check(nf90_put_var(ncid, var_1_id, cons))
    elseif (outtype == 'coor') then 
        call check(nf90_put_var(ncid, var_1_id, x))
    elseif (outtype == 'ien') then 
        call check(nf90_put_var(ncid, var_1_id, ien))
    elseif (outtype == 'nsmp') then 
        call check(nf90_put_var(ncid, var_1_id, nsmp_tmp))
    endif 
    ! Close the file.
    call check(nf90_close(ncid))
    !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

end subroutine netcdf_write

! #2
! netcdf_write_on_fault writes on-fault quantities to netcdf files.
!
! Each fault has its own local (strike, dip) grid, sized from its own
! fltxyz(:,:,ift) box, not from the whole-model mesh (nxt, nzt): those only
! happen to coincide with the fault's own extent when a single fault spans
! the entire model, which is the only case exercised so far. A fault node's
! local (i,j) is recovered from its own coordinates rather than assumed from
! loop order, since node ordinal number is not guaranteed to enumerate a
! fault's own strike-major grid once a fault does not span the whole mesh.
! The written arrays are padded to the largest fault's extent
! (nfxMax, nfzMax) and carry an explicit nid_fault dimension of size ntotft,
! so the same statements run whether ntotft is 1 or more.
subroutine netcdf_write_on_fault(outfile)
    use netcdf
    use globalvar
    implicit none
    character (len = 50 ) :: outfile, lat_name, lon_name, flt_name, lat_units, lon_units, flt_units, UNITS
    character (len = 50), allocatable, dimension(:) :: var_name, var_unit
    integer (kind = 4) :: ncid, lat_dimid, lon_dimid, flt_dimid, lat_varid, lon_varid, flt_varid, &
        var_id(20), i, j, n, ift, nlat, nlon, nvar
    integer (kind = 4) :: dimids(3)
    integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index, flt_index
    integer (kind = 4), allocatable, dimension(:) :: nfx, nfz
    integer (kind = 4) :: nfxMax, nfzMax
    real (kind = dp), allocatable, dimension(:,:,:,:) :: on_fault_vars

    nvar = 15 ! nvar variables.

    ! Per-fault local grid extent, from that fault's own bounding box.
    allocate(nfx(ntotft), nfz(ntotft))
    do ift = 1, ntotft
        nfx(ift) = int((fltxyz(2,1,ift) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
        nfz(ift) = int((fltxyz(2,3,ift) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
    enddo
    nfxMax = maxval(nfx)
    nfzMax = maxval(nfz)

    allocate(var_name(nvar), var_unit(nvar), on_fault_vars(nfxMax, nfzMax, nvar, ntotft)) ! on_fault_vars(lon,lat,nvar,fault)
    on_fault_vars = 0.0d0

    UNITS     = 'units'
    lat_name  = 'nid_dip'
    lon_name  = 'nid_strike'
    flt_name  = 'nid_fault'
    lat_units = 'unit'
    lon_units = 'unit'
    flt_units = 'unit'
    var_name  = [ character(len=20) :: 'shear_strike', 'shear_dip', 'effective_normal', 'slip_rate' , 'state_variable', &
        'state_normal', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs', 'slips', 'slipd', 'slipn']
    var_unit  = [ character(len=20) :: 'Pa'          , 'Pa'       , 'Pa'              , 'm/s'       , 'unit'          , &
        'Pa'          , 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'm/s', 'slips', 'slipd', 'slipn']
    nlat = nfzMax ! Total nodes along dip.
    nlon = nfxMax ! Total nodes along strike.

    allocate(lat_index(nlat), lon_index(nlon), flt_index(ntotft))

    lat_index = (/ (i, i = 1, nlat) /)
    lon_index = (/ (i, i = 1, nlon) /)
    flt_index = (/ (i, i = 1, ntotft) /)

    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            i = int((x(1,nsmp(1,n,ift)) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
            j = int((x(3,nsmp(1,n,ift)) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
            on_fault_vars(i,j,1,ift)  = fric(28, n, ift) ! tstk0
            on_fault_vars(i,j,2,ift)  = fric(29, n, ift) ! tdip0
            on_fault_vars(i,j,3,ift)  = fric(30, n, ift) ! tnorm0
            on_fault_vars(i,j,4,ift)  = fric(26, n, ift) ! sliprate
            on_fault_vars(i,j,5,ift)  = fric(20, n, ift) ! state
            on_fault_vars(i,j,6,ift)  = fric(23, n, ift) ! state variable for normal stress, theta_pc
            on_fault_vars(i,j,7,ift)  = fric(31, n, ift) ! vxm
            on_fault_vars(i,j,8,ift)  = fric(32, n, ift) ! vym
            on_fault_vars(i,j,9,ift)  = fric(33, n, ift) ! vzm
            on_fault_vars(i,j,10,ift) = fric(34, n, ift) ! vxs
            on_fault_vars(i,j,11,ift) = fric(35, n, ift) ! vys
            on_fault_vars(i,j,12,ift) = fric(36, n, ift) ! vzs
            on_fault_vars(i,j,13,ift) = fric(71, n, ift) ! vxs
            on_fault_vars(i,j,14,ift) = fric(72, n, ift) ! vys
            on_fault_vars(i,j,15,ift) = fric(73, n, ift) ! vzs
        enddo
    enddo
    ! Create the netCDF file.
    call check(nf90_create(outfile, NF90_CLOBBER, ncid))

    ! Define the dimensions.
    call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
    call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
    call check(nf90_def_dim(ncid, flt_name, ntotft, flt_dimid))

    ! Define coordiante variables. They will hold the coordinate
    ! information, that is, the latitudes (y), and longitudes (x). A varid is
    ! returned for each.
    call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
    call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
    call check(nf90_def_var(ncid, flt_name, NF90_INT, flt_dimid, flt_varid))

    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
    call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )
    call check( nf90_put_att(ncid, flt_varid, UNITS, flt_units) )

    ! Define the netcdf variables. The dimids array is used to pass the
    ! dimids of the dimensions of the netCDF variables.
    dimids = (/ lon_dimid, lat_dimid, flt_dimid /)

    do i = 1, nvar
        call check(nf90_def_var(ncid, var_name(i), NF90_REAL, dimids, var_id(i)))
        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        call check(nf90_put_att(ncid, var_id(i), UNITS, var_unit(i)))
    enddo

    ! End definitions.
    call check(nf90_enddef(ncid))

    ! Write data.
    ! Write the coordinate variable data. This will put the x, and y
    ! of our data grid into the netCDF file.
    call check(nf90_put_var(ncid, lat_varid, lat_index))
    call check(nf90_put_var(ncid, lon_varid, lon_index))
    call check(nf90_put_var(ncid, flt_varid, flt_index))

    ! Write the data. This will write our displacement fields "cons" which is defined
    ! globally. Its dimension is 3 by numnp (row by column).
    do i = 1,nvar
        call check(nf90_put_var(ncid, var_id(i), on_fault_vars(:,:,i,:)))
    enddo

    ! Close the file.
    call check(nf90_close(ncid))
    !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

end subroutine netcdf_write_on_fault

! #3
subroutine netcdf_write_roughness(outfile)
    use netcdf
    use globalvar
    implicit none 
    character (len = 50 ) :: outfile, outtype, lat_name, lon_name, lat_units, lon_units, UNITS
    character (len = 50), allocatable, dimension(:) :: var_name, var_unit
    integer (kind = 4) :: ncid, lat_dimid, lon_dimid, lat_varid, lon_varid, var_id(20), ilat, ilon, i, j, nlat, nlon, nvar
    integer (kind = 4) :: dimids(2)
    integer (kind = 4), allocatable, dimension(:) :: lat_index, lon_index
    real (kind = dp), allocatable, dimension(:,:,:) :: on_fault_vars

    nvar = 3 ! nvar variables.
    allocate(var_name(nvar), var_unit(nvar), on_fault_vars(nnx, nnz, nvar)) ! on_fault_vars(lon,lat,nvar)
    on_fault_vars = 0.0d0
    
    UNITS     = 'units'
    lat_name  = 'nid_dip'
    lon_name  = 'nid_strike'
    lat_units = 'unit'
    lon_units = 'unit'
    var_name  = [ character(len=20) :: 'peak', 'pypx', 'pypz' ]
    var_unit  = [ character(len=20) :: 'm'          , 'unit'       , 'unit']
    nlat      = nnz ! Total nodes along dip. 
    nlon      = nnx ! Total nodes along strike.  
    
    allocate(lat_index(nlat))
    allocate(lon_index(nlon))
    
    lat_index = (/ (i, i = 1, nlat) /)
    lon_index = (/ (i, i = 1, nlon) /)
    
    do i = 1, nnx
        do j = 1, nnz
            on_fault_vars(i,j,1) = rough_geo(1, (i-1)*nnz+j) ! tstk0
            on_fault_vars(i,j,2) = rough_geo(2, (i-1)*nnz+j) ! tdip0
            on_fault_vars(i,j,3) = rough_geo(3, (i-1)*nnz+j) ! tnorm0
        enddo 
    enddo 
    ! Create the netCDF file.
    call check(nf90_create(outfile, NF90_CLOBBER, ncid))
    
    ! Define the dimensions.
    call check(nf90_def_dim(ncid, lat_name, nlat, lat_dimid))
    call check(nf90_def_dim(ncid, lon_name, nlon, lon_dimid))
    
    ! Define coordiante variables. They will hold the coordinate 
    ! information, that is, the latitudes (y), and longitudes (x). A varid is 
    ! returned for each.
    call check(nf90_def_var(ncid, lat_name, NF90_INT, lat_dimid, lat_varid))
    call check(nf90_def_var(ncid, lon_name, NF90_INT, lon_dimid, lon_varid))
    
    ! Assign units attributes to coordinate var data. This attaches a
    ! text attribute to each of the coordinate variables, containing the
    ! units.
    call check( nf90_put_att(ncid, lat_varid, UNITS, lat_units) )
    call check( nf90_put_att(ncid, lon_varid, UNITS, lon_units) )    

    ! Define the netcdf variables. The dimids array is used to pass the 
    ! dimids of the dimensions of the netCDF variables.
    dimids = (/ lon_dimid, lat_dimid /)
    
    do i = 1, nvar
        call check(nf90_def_var(ncid, var_name(i), NF90_REAL, dimids, var_id(i))) 
        ! Assign units attributes to the pressure and temperature netCDF
        ! variables.
        call check(nf90_put_att(ncid, var_id(i), UNITS, var_unit(i)))
    enddo 

    ! End definitions.
    call check(nf90_enddef(ncid))
    
    ! Write data.
    ! Write the coordinate variable data. This will put the x, and y 
    ! of our data grid into the netCDF file.
    call check(nf90_put_var(ncid, lat_varid, lat_index))
    call check(nf90_put_var(ncid, lon_varid, lon_index))
    
    ! Write the data. This will write our displacement fields "cons" which is defined 
    ! globally. Its dimension is 3 by numnp (row by column).
    do i = 1,nvar
        call check(nf90_put_var(ncid, var_id(i), on_fault_vars(:,:,i)))
    enddo
    
    ! Close the file.
    call check(nf90_close(ncid))
    !call writegrid(outfile, xpos, ypos, zpos, data_arr, nxtmp, nytmp, nztmp)

end subroutine netcdf_write_roughness

! Subroutine #4.
! netcdf_read_on_fault reads in on-fault quantities from netcdf files created by case.setup.
!
! The file may be in either of two layouts, detected from the rank of its
! variables -- never from ntotft:
!   - legacy 2-D layout (dims: dip, strike), written by every compset before
!     multi-fault support, for the one fault that then spanned the whole
!     model grid.
!   - 3-D layout (dims: dip, strike, nid_fault), written by case.setup's
!     multi-fault netcdf_write_on_fault_vars, padded to the largest fault's
!     (nfxMax, nfzMax) extent -- the same convention netcdf_write_on_fault
!     (this file) uses for its own output, see that subroutine's header.
! Either way, a fault node's local (i,j) is recovered from its own
! coordinates rather than assumed from loop order (node ordinal number is
! not guaranteed to enumerate a fault's own strike-major grid once a fault
! does not span the whole mesh), so this same code runs whether ntotft is 1
! or more.
subroutine netcdf_read_on_fault(infile)
    use netcdf
    use globalvar
    implicit none
    character (len = 50 ) :: infile
    integer (kind = 4) :: ncid, var_id(20), i, j, n, ift, nvar, ndims_var
    integer (kind = 4) :: dimids3(3), dimlen1, dimlen2, dimlen3
    integer (kind = 4), allocatable, dimension(:) :: nfx, nfz
    integer (kind = 4) :: nfxMax, nfzMax
    real (kind = dp), allocatable, dimension(:,:,:,:) :: on_fault_vars

    nvar = 9

    ! Per-fault local grid extent, from that fault's own bounding box; same
    ! formula as netcdf_write_on_fault, so the padded shape agrees exactly.
    allocate(nfx(ntotft), nfz(ntotft))
    do ift = 1, ntotft
        nfx(ift) = int((fltxyz(2,1,ift) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
        nfz(ift) = int((fltxyz(2,3,ift) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
    enddo
    nfxMax = maxval(nfx)
    nfzMax = maxval(nfz)

    allocate(on_fault_vars(nfxMax,nfzMax,ntotft,nvar))
    on_fault_vars = 0.0d0

    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    call check( nf90_open(infile, NF90_NOWRITE, ncid))

    ! Get the varid of the data variables, based on their names.
    call check( nf90_inq_varid(ncid, "a",  var_id(1)))
    call check( nf90_inq_varid(ncid, "b",  var_id(2)))
    call check( nf90_inq_varid(ncid, "Dc", var_id(3)))
    call check( nf90_inq_varid(ncid, "v0", var_id(4)))
    call check( nf90_inq_varid(ncid, "r0", var_id(5)))
    call check( nf90_inq_varid(ncid, "init_slip_rate",     var_id(6)))
    call check( nf90_inq_varid(ncid, "init_shear_stress",  var_id(7)))
    call check( nf90_inq_varid(ncid, "init_normal_stress", var_id(8)))
    call check( nf90_inq_varid(ncid, "init_state",         var_id(9)))

    ! Detect the file layout from the rank of "a": 2 for the legacy
    ! single-fault format, 3 once a nid_fault dimension is present.
    dimids3 = 0
    call check( nf90_inquire_variable(ncid, var_id(1), ndims = ndims_var, dimids = dimids3))

    if (ndims_var == 2) then
        ! Legacy format has no fault dimension at all -- it structurally
        ! cannot supply more than fault 1's data. Fail loudly rather than
        ! silently leaving fault 2..ntotft zeroed.
        if (ntotft > 1) then
            print *, 'netcdf_read_on_fault: ', trim(infile), &
                ' has no nid_fault dimension (rank-2 legacy format) but ntotft = ', &
                ntotft, '. Regenerate on_fault_vars_input.nc with case.setup.'
            stop 'Stopped'
        endif
        call check( nf90_inquire_dimension(ncid, dimids3(1), len = dimlen1))
        call check( nf90_inquire_dimension(ncid, dimids3(2), len = dimlen2))
        if (dimlen1 /= nfx(1) .or. dimlen2 /= nfz(1)) then
            print *, 'netcdf_read_on_fault: ', trim(infile), &
                ' strike/dip extent (', dimlen1, dimlen2, &
                ') does not match fault 1''s mesh extent (', nfx(1), nfz(1), ').'
            stop 'Stopped'
        endif
        do i = 1, nvar
            call check( nf90_get_var(ncid, var_id(i), on_fault_vars(1:nfx(1), 1:nfz(1), 1, i)))
        enddo
    elseif (ndims_var == 3) then
        call check( nf90_inquire_dimension(ncid, dimids3(1), len = dimlen1))
        call check( nf90_inquire_dimension(ncid, dimids3(2), len = dimlen2))
        call check( nf90_inquire_dimension(ncid, dimids3(3), len = dimlen3))
        if (dimlen1 /= nfxMax .or. dimlen2 /= nfzMax .or. dimlen3 /= ntotft) then
            print *, 'netcdf_read_on_fault: ', trim(infile), &
                ' shape (', dimlen1, dimlen2, dimlen3, &
                ') does not match expected (nfxMax, nfzMax, ntotft) = (', &
                nfxMax, nfzMax, ntotft, ').'
            stop 'Stopped'
        endif
        do i = 1, nvar
            call check( nf90_get_var(ncid, var_id(i), on_fault_vars(1:nfxMax, 1:nfzMax, 1:ntotft, i)))
        enddo
    else
        print *, 'netcdf_read_on_fault: ', trim(infile), &
            ' variable "a" has unsupported rank ', ndims_var, '; expected 2 or 3.'
        stop 'Stopped'
    endif

    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            i = int((x(1,nsmp(1,n,ift)) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
            j = int((x(3,nsmp(1,n,ift)) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
            fric(9,  n, ift) = on_fault_vars(i,j,ift,1) ! a
            fric(10, n, ift) = on_fault_vars(i,j,ift,2) ! b
            fric(11, n, ift) = on_fault_vars(i,j,ift,3) ! Dc
            fric(12, n, ift) = on_fault_vars(i,j,ift,4) ! v0
            fric(13, n, ift) = on_fault_vars(i,j,ift,5) ! r0
            fric(46, n, ift) = on_fault_vars(i,j,ift,6) ! init_slip_rate
            fric(8,  n, ift) = on_fault_vars(i,j,ift,7) ! shear
            fric(7,  n, ift) = on_fault_vars(i,j,ift,8) ! norm
            fric(20, n, ift) = on_fault_vars(i,j,ift,9) ! state variable
            ! The aging law (friclaw 3) carries theta, a time in seconds. The
            ! slip law (friclaw 4) carries psi = f* + b*ln(V* theta / Dc), which
            ! is dimensionless and of order f*. Compsets specify the initial
            ! state as theta, following SEAS BP8 eq. (30), so convert it here
            ! rather than handing the slip law a number 10^9 times too large:
            ! exp(psi/a) then overflows and the Newton solve returns NaN on the
            ! very first step.
            if (friclaw == 4) then
                fric(20, n, ift) = fric(13, n, ift) &
                    + fric(10, n, ift) * dlog(fric(12, n, ift) &
                    * on_fault_vars(i,j,ift,9) / fric(11, n, ift))
            endif
            ! Snapshot the initial state. The t = 0 row of the station files is
            ! written when the file is created, which happens at the end of the
            ! run, so reading fric(20) there reports the *final* state under the
            ! label t = 0. Keep the initial value so that row is honest.
            fric(48, n, ift) = fric(20, n, ift) ! initial state
            fric(47, n, ift) = fric(46, n, ift) ! peak slip rate
            fric(23, n, ift) = abs(fric(7, n, ift)) ! initialize theta_pc as abs(normal stress)
        enddo
    enddo

    ! Close the file, freeing all resources.
    call check( nf90_close(ncid))

end subroutine netcdf_read_on_fault

! Subroutine #5.
! netcdf_read_on_fault_restart reads in additional on-fault quantities from netcdf files created by previous cycles.
!
! Phase one re-reads a, b, Dc, v0, r0 from infile1 (on_fault_vars_input.nc,
! the same file netcdf_read_on_fault reads for cycle 1): rank-detected and
! coordinate-mapped exactly as there, see that subroutine's header. Phase
! two reads the previous cycle's state from infile2 (fault.r.nc, written by
! netcdf_write_on_fault) and was already ift-looped.
subroutine netcdf_read_on_fault_restart(infile1, infile2)
    use netcdf
    use globalvar
    implicit none
    character (len = 50 ) :: infile1, infile2
    integer (kind = 4)    :: ncid,  var_id(20), i, j, n, ift, nvar, ndims_var
    integer (kind = 4) :: dimids3(3), dimlen1, dimlen2, dimlen3
    integer (kind = 4), allocatable, dimension(:) :: nfx, nfz
    integer (kind = 4) :: nfxMax, nfzMax
    real (kind = dp), allocatable, dimension(:,:,:,:) :: on_fault_vars
    real (kind = dp), allocatable, dimension(:,:,:,:) :: on_fault_vars4

    ! Per-fault local grid extent, shared by both phases below.
    allocate(nfx(ntotft), nfz(ntotft))
    do ift = 1, ntotft
        nfx(ift) = int((fltxyz(2,1,ift) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
        nfz(ift) = int((fltxyz(2,3,ift) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
    enddo
    nfxMax = maxval(nfx)
    nfzMax = maxval(nfz)

    ! Read in 5 variables a, b, Dc, v0, r0.
    nvar = 5
    allocate(on_fault_vars(nfxMax,nfzMax,ntotft,nvar))
    on_fault_vars = 0.0d0

    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    call check( nf90_open(infile1, NF90_NOWRITE, ncid))

    ! Get the varid of the data variables, based on their names.
    call check( nf90_inq_varid(ncid, "a",  var_id(1)))
    call check( nf90_inq_varid(ncid, "b",  var_id(2)))
    call check( nf90_inq_varid(ncid, "Dc", var_id(3)))
    call check( nf90_inq_varid(ncid, "v0", var_id(4)))
    call check( nf90_inq_varid(ncid, "r0", var_id(5)))

    ! Detect the file layout from the rank of "a", as netcdf_read_on_fault does.
    dimids3 = 0
    call check( nf90_inquire_variable(ncid, var_id(1), ndims = ndims_var, dimids = dimids3))

    if (ndims_var == 2) then
        if (ntotft > 1) then
            print *, 'netcdf_read_on_fault_restart: ', trim(infile1), &
                ' has no nid_fault dimension (rank-2 legacy format) but ntotft = ', &
                ntotft, '. Regenerate on_fault_vars_input.nc with case.setup.'
            stop 'Stopped'
        endif
        call check( nf90_inquire_dimension(ncid, dimids3(1), len = dimlen1))
        call check( nf90_inquire_dimension(ncid, dimids3(2), len = dimlen2))
        if (dimlen1 /= nfx(1) .or. dimlen2 /= nfz(1)) then
            print *, 'netcdf_read_on_fault_restart: ', trim(infile1), &
                ' strike/dip extent (', dimlen1, dimlen2, &
                ') does not match fault 1''s mesh extent (', nfx(1), nfz(1), ').'
            stop 'Stopped'
        endif
        do i = 1, nvar
            call check( nf90_get_var(ncid, var_id(i), on_fault_vars(1:nfx(1), 1:nfz(1), 1, i)))
        enddo
    elseif (ndims_var == 3) then
        call check( nf90_inquire_dimension(ncid, dimids3(1), len = dimlen1))
        call check( nf90_inquire_dimension(ncid, dimids3(2), len = dimlen2))
        call check( nf90_inquire_dimension(ncid, dimids3(3), len = dimlen3))
        if (dimlen1 /= nfxMax .or. dimlen2 /= nfzMax .or. dimlen3 /= ntotft) then
            print *, 'netcdf_read_on_fault_restart: ', trim(infile1), &
                ' shape (', dimlen1, dimlen2, dimlen3, &
                ') does not match expected (nfxMax, nfzMax, ntotft) = (', &
                nfxMax, nfzMax, ntotft, ').'
            stop 'Stopped'
        endif
        do i = 1, nvar
            call check( nf90_get_var(ncid, var_id(i), on_fault_vars(1:nfxMax, 1:nfzMax, 1:ntotft, i)))
        enddo
    else
        print *, 'netcdf_read_on_fault_restart: ', trim(infile1), &
            ' variable "a" has unsupported rank ', ndims_var, '; expected 2 or 3.'
        stop 'Stopped'
    endif

    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            i = int((x(1,nsmp(1,n,ift)) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
            j = int((x(3,nsmp(1,n,ift)) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
            fric(9,  n, ift) = on_fault_vars(i,j,ift,1) ! a
            fric(10, n, ift) = on_fault_vars(i,j,ift,2) ! b
            fric(11, n, ift) = on_fault_vars(i,j,ift,3) ! Dc
            fric(12, n, ift) = on_fault_vars(i,j,ift,4) ! v0
            fric(13, n, ift) = on_fault_vars(i,j,ift,5) ! r0
        enddo
    enddo
    ! Close the file, freeing all resources.
    call check( nf90_close(ncid))

    deallocate(on_fault_vars)

    ! Phase two, read in initial conditions from restart files fault.r.nc
    ! (written by netcdf_write_on_fault, which carries an explicit nid_fault
    ! dimension of size ntotft and per-fault local (strike, dip) grids -- see
    ! that subroutine's header comment).
    nvar = 12
    allocate(on_fault_vars4(nfxMax,nfzMax,nvar,ntotft))
    ! Open the file. NF90_NOWRITE tells netCDF we want read-only access to the file.
    call check( nf90_open(infile2, NF90_NOWRITE, ncid))

    ! Get the varid of the data variables, based on their names.
    ! 'shear_strike', 'shear_dip', 'effective_normal', 'slip_rate' , 'state_variable', 'vxm', 'vym', 'vzm', 'vxs', 'vys', 'vzs'
    call check( nf90_inq_varid(ncid, "shear_strike",     var_id(1)))
    call check( nf90_inq_varid(ncid, "shear_dip",        var_id(2)))
    call check( nf90_inq_varid(ncid, "effective_normal", var_id(3)))
    call check( nf90_inq_varid(ncid, "slip_rate",        var_id(4)))
    call check( nf90_inq_varid(ncid, "state_variable",   var_id(5)))
    call check( nf90_inq_varid(ncid, "state_normal",     var_id(6)))
    call check( nf90_inq_varid(ncid, "vxm", var_id(7) ))
    call check( nf90_inq_varid(ncid, "vym", var_id(8) ))
    call check( nf90_inq_varid(ncid, "vzm", var_id(9) ))
    call check( nf90_inq_varid(ncid, "vxs", var_id(10)))
    call check( nf90_inq_varid(ncid, "vys", var_id(11)))
    call check( nf90_inq_varid(ncid, "vzs", var_id(12)))
    ! Read the data
    do i = 1, nvar
        call check( nf90_get_var(ncid, var_id(i), on_fault_vars4(:,:,i,:)))
    enddo

    do ift = 1, ntotft
        do n = 1, nftnd(ift)
            i = int((x(1,nsmp(1,n,ift)) - fltxyz(1,1,ift))/dx + 0.5d0) + 1
            j = int((x(3,nsmp(1,n,ift)) - fltxyz(1,3,ift))/dx + 0.5d0) + 1
            fric(8,  n, ift) = on_fault_vars4(i,j,1,ift)! tstk0
            fric(49, n, ift) = on_fault_vars4(i,j,2,ift)! tdip0
            fric(7,  n, ift) = on_fault_vars4(i,j,3,ift)! tnorm0
            fric(46, n, ift) = on_fault_vars4(i,j,4,ift)! sliprate
            fric(20, n, ift) = on_fault_vars4(i,j,5,ift)! state
            fric(23, n, ift) = on_fault_vars4(i,j,6,ift)! state_normal
            fric(31, n, ift) = on_fault_vars4(i,j,7,ift)! vxm
            fric(32, n, ift) = on_fault_vars4(i,j,8,ift)! vym
            fric(33, n, ift) = on_fault_vars4(i,j,9,ift)! vzm
            fric(34, n, ift) = on_fault_vars4(i,j,10,ift)! vxs
            fric(35, n, ift) = on_fault_vars4(i,j,11,ift)! vys
            fric(36, n, ift) = on_fault_vars4(i,j,12,ift)! vzs
            fric(48, n, ift) = fric(20, n, ift)! state at cycle start
            fric(47, n, ift) = fric(46, n, ift)! peak slip rate
            !fric(23, n, ift) = abs(fric(7, n, ift))! initialize theta_pc as abs(normal stress)
        enddo
    enddo
    ! Close the file, freeing all resources.
    call check( nf90_close(ncid))

    deallocate(on_fault_vars4)

end subroutine netcdf_read_on_fault_restart

subroutine check(status)
    use netcdf
    integer, intent ( in) :: status
    
    if(status /= nf90_noerr) then 
      print *, nf90_strerror(status)
      stop "Stopped"
    end if
end subroutine check  