subroutine output_onfault_st

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j, k

    if(n4onf > 0) then
        do i = 1,n4onf
            j = anonfs(3,i)
            if (j==1) then  !main fault stations
                sttmp = '      '
                dptmp = '      '
                if (bp == 8) then
                    ! SEAS BP8 names stations in metres with an explicit sign,
                    ! e.g. fltst_strk+000dp-200.
                    write(sttmp,'(SP,i4.3)') int(xonfs(1,anonfs(2,i),j))
                    write(dptmp,'(SP,i4.3)') int(-xonfs(2,anonfs(2,i),j))
                elseif (bp == 7) then
                    write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j))
                    write(dptmp,'(i4.3)') int(-xonfs(2,anonfs(2,i),j))
                else
                    write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)/1000.d0)
                    write(dptmp,'(i4.3)') int((-xonfs(2,anonfs(2,i),j))/1000.d0)
                endif
                open(51,file='fltst_strk'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')     
            endif
            !write(51,*) '# This is the file header'
            !write(51,*) '# problem = San-Ti'
            !write(51,*) '# Author = Sophon'
            !write(51,*) '# code = EQquasi'
            !write(51,*) '# version = 1.x'   
            if (bp == 8) then
                ! SEAS BP8 time series: 11 fields, see section 4.1 of the
                ! benchmark description.
                write(51,'(A)') '# problem=SEAS Benchmark BP8-QD-'//trim(adjustl(bp8tag))
                write(51,'(A)') '# code=EQquasi'
                write(51,'(A)') '# version='//EQQUASI_VERSION
                write(51,'(A)') '# modeler=D.Liu'
                write(51,'(A)') '# date='//trim(adjustl(runDate))
                write(51,'(A,E15.7,A)') '# element_size=', dx, ' m'
                write(51,'(A)') '# location= on fault'
                write(51,'(A)') '# Column #1 = Time (s)'
                write(51,'(A)') '# Column #2 = Slip_2 (m)'
                write(51,'(A)') '# Column #3 = Slip_3 (m)'
                write(51,'(A)') '# Column #4 = Slip_rate_2 (log10 m/s)'
                write(51,'(A)') '# Column #5 = Slip_rate_3 (log10 m/s)'
                write(51,'(A)') '# Column #6 = Shear_stress_2 (MPa)'
                write(51,'(A)') '# Column #7 = Shear_stress_3 (MPa)'
                write(51,'(A)') '# Column #8 = Pore_pressure (MPa)'
                write(51,'(A)') '# Column #9 = Darcy_velocity_2 (m/s)'
                write(51,'(A)') '# Column #10 = Darcy_velocity_3 (m/s)'
                write(51,'(A)') '# Column #11 = State (log10 s)'
                write(51,'(A)') '# The line below lists the names of the data fields'
                write(51,'(A)') 't slip_2 slip_3 slip_rate_2 slip_rate_3 shear_stress_2 '// &
                                'shear_stress_3 pore_pressure darcy_vel_2 darcy_vel_3 state'
                write(51,'(A)') '# Here is the time-series data.'
                ! Report the prescribed initial condition at t = 0. fltsta only
                ! holds computed steps, so the series would otherwise begin at
                ! the first step and carry no t = 0 record at all.
                k = anonfs(1,i)
                write(51,'(E22.14,10E15.7)') 0.0d0,        & ! t
                    0.0d0, 0.0d0,                          & ! slip_2, slip_3
                    dlog10(max(fric(46,k,1), 1.0d-30)),    & ! log10 V_init
                    dlog10(1.0d-20),                       & ! log10 V_zero
                    fric(8,k,1)/1.d6, 0.0d0,               & ! tau^0, tau_3
                    0.0d0, 0.0d0, 0.0d0,                   & ! p, q_2, q_3
                    dlog10(max(fric(11,k,1)/max(fric(46,k,1),1.0d-30), 1.0d-30))
                do j = 1,it-1
                    write(51,'(E22.14,10E15.7)') fltsta(1,j,i),   & ! Time in sec
                        fltsta(5,j,i),                            & ! slip_2, along strike
                        -fltsta(6,j,i),                           & ! slip_3, positive downwards
                        dlog10(max(abs(fltsta(14,j,i)), 1.0d-30)),& ! log10 |slip rate along strike|
                        dlog10(max(abs(fltsta(3,j,i)),  1.0d-30)),& ! log10 |slip rate along dip|
                        fltsta(8,j,i)/1.d6,                       & ! shear stress along strike, MPa
                        -fltsta(9,j,i)/1.d6,                      & ! shear stress along dip, MPa
                        fltsta(11,j,i)/1.d6,                      & ! pore pressure change, MPa
                        fltsta(12,j,i),                           & ! Darcy velocity along strike, m/s
                        -fltsta(13,j,i),                          & ! Darcy velocity along dip, m/s
                        dlog10(max(fltsta(4,j,i), 1.0d-30))         ! log10 state variable
                enddo
                close(51)
                cycle
            endif

            do j = 1,it-1
                write(51,'( 9e32.21e4)') fltsta(1,j,i), & ! Time in sec
                    fltsta(5,j,i), & ! slip_strike, along x;
                    -fltsta(6,j,i), & ! -slip_dip, along -z, downwards;
                    fltsta(2,j,i), & ! trial slip rate magnitude, m/s;
                    fltsta(3,j,i), & ! slip rate along dip
                    fltsta(8,j,i)/1.d6, & ! shear stress along strike, MPa
                    -fltsta(9,j,i)/1.d6, & ! minus shear stress along dip, MPa
                    fltsta(10,j,i)/1.d6, & ! effective normal stress, MPa; negative compressive. 
                    dlog10(fltsta(4,j,i)) ! log 10 of state variable. 
                enddo
            close(51)
        enddo
    endif
end subroutine output_onfault_st

subroutine output_offfault_st
    
    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j     
    character (len = 30) :: bodytmp
    
    if(n4out>0) then
        do i=1,n4out
            bodytmp = '      '
            sttmp = '      '
            dptmp = '      '
            if (bp == 7 .or. bp == 8) then
                ! Metre-scale models: dividing by 1000 collapses distinct
                ! stations onto the same filename and silently overwrites them.
                write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i)))
                write(sttmp,'(i4.3)')   int(x4nds(1,an4nds(1,i)))
                write(dptmp,'(i4.3)')   int(abs(x4nds(3,an4nds(1,i))))
            else
                write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i))/1000.d0)
                write(sttmp,'(i4.3)')   int(x4nds(1,an4nds(1,i))/1000.d0)
                write(dptmp,'(i4.3)')   int(x4nds(3,an4nds(1,i))/1000.d0)
            endif
            open(51,file='srfst_strk'//trim(adjustl(sttmp))//'st'//trim(adjustl(bodytmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

            !write(51,*) '# This is the file header'
            !write(51,*) '# problem = San-Ti'
            !write(51,*) '# Author = Sophon'
            !write(51,*) '# code = EQquasi'
            !write(51,*) '# version = 1.x'   
            !write(51,*) '# Time series in 11 columns in format E15.7 for data'
            !write(51,*) '# Column #1 = Time (s)'
            !write(51,*) '# Column #2 = horizontal displacement (m)'
            !write(51,*) '# Column #3 = horizontal velocity (m/s)'
            !write(51,*) '# Column #4 = vertical displacement (m)'
            !write(51,*) '# Column #5 = vertical velocity (m/s)'
            !write(51,*) '# Column #6 = normal displacement (m)'
            !write(51,*) '# Column #7 = normal velocity (m/s)'
            !write(51,*) '#'
            !write(51,*) '# The line below lists the names of the data fields:'
            !write(51,*) 't h-disp h-vel v-disp v-vel n-disp n-vel'
            do j=1,it-1
                write(51,'( E21.13,6E15.7)') dout(1,j),dout((i-1)*6+2,j), &
                dout((i-1)*6+3,j),-dout((i-1)*6+6,j),-dout((i-1)*6+7,j), &
                dout((i-1)*6+4,j),dout((i-1)*6+5,j)
            enddo
            close(51)
        enddo
    endif
end subroutine output_offfault_st

subroutine output_onfault_transfer

    use globalvar
    implicit none
    
    integer (kind = 4) :: i
    
    if(nftnd(1) > 0) then
        open(unit=1111,file='cplot_EQquasi.txt',status='unknown')    !rupture time
        ! write(1111,*) '# This is the file header:'
        ! write(1111,*) '# problem = San-Ti'
        ! write(1111,*) '# author = sophon'
        ! write(1111,*) '# date = 2424/01/01'
        ! write(1111,*) '# code = Qquasi'
        ! write(1111,*) '# version = 1.x'
        write(1111,'(1x,16e32.21e4)') (x(1,nsmp(1,i,1)), -x(3,nsmp(1,i,1)), & ! xcoor, -zcoor
            fric(26,i,1), & ! trial slip rate magnitude, m/s;
            fric(20,i,1), & ! state variable in RSF;
            fric(28,i,1), & ! shear stress along strike, Pa;
            fric(29,i,1), & ! shear stress along dip, Pa;
            fric(30,i,1), & ! effective normal stress, Pa;
            fric(31,i,1), fric(32,i,1), fric(33,i,1), & ! master node nodal velocity along x, y, z; master node on y+ side of the fault. 
            fric(34,i,1), fric(35,i,1), fric(36,i,1), & ! slave node nodal velocity along x, y, z; slave node on the y- side of the fault. 
            fric(44,i,1), fric(45,i,1), & ! empty for now. 
            fnft(i,1), & ! rupture time at this loc.
                i = 1,nftnd(1))    
        close(1111)
    endif
    
end subroutine output_onfault_transfer

subroutine output_timedy
    use globalvar
    implicit none
    
    open(1112,file='tdyna.txt',form = 'formatted', status = 'unknown')
        write(1112,'(2e32.21e4)') tdynastart, tdynaend
    close(1112)

end subroutine output_timedy

subroutine output_globaldat
    use globalvar
    implicit none
    integer (kind = 4) :: i 
    
    if (bp == 8) then
        ! SEAS BP8 source parameter time series, see section 4.2.
        open(1113,file='global.dat',form='formatted',status='unknown')
            write(1113,'(A)') '# problem=SEAS Benchmark BP8-QD-'//trim(adjustl(bp8tag))
            write(1113,'(A)') '# code=EQquasi'
            write(1113,'(A)') '# version='//EQQUASI_VERSION
            write(1113,'(A)') '# modeler=D.Liu'
            write(1113,'(A)') '# date='//trim(adjustl(runDate))
            write(1113,'(A,E15.7,A)') '# element_size=', dx, ' m'
            write(1113,'(A)') '# location= frictional domain'
            write(1113,'(A)') '# Column #1 = Time (s)'
            write(1113,'(A)') '# Column #2 = Max_slip_rate (log10 m/s)'
            write(1113,'(A)') '# Column #3 = Moment_density_rate (N/s)'
            write(1113,'(A)') '# The line below lists the names of the data fields'
            write(1113,'(A)') 't max_slip_rate moment_density_rate'
            write(1113,'(A)') '# Here is the time-series data.'
            ! Initial condition at t = 0, as for the station files.
            write(1113,'(E22.14,2E15.7)') 0.0d0, &
                dlog10(max(maxval(fric(46,1:nftnd(1),1)), 1.0d-30)), 0.0d0
            do i = 1, it-1
                write(1113,'(E22.14,2E15.7)') globaldat(1,i), &
                    dlog10(max(globaldat(2,i), 1.0d-30)), &
                    globaldat(3,i)
            enddo
        close(1113)
        return
    endif

    open(1113,file='global.dat',form='formatted',status='unknown')
        write(1113,'(7e32.21e4)') (globaldat(1,i), &
            globaldat(2,i), &
            globaldat(3,i), &
            globaldat(4,i), &
            globaldat(5,i), &
            globaldat(6,i), &
            globaldat(7,i), i = 1,it-1)
    close(1113)

end subroutine output_globaldat


subroutine output_prof

    use globalvar
    implicit none
    
    integer (kind = 4)  :: i, ift 
    
    ift = 1 
    if (bp == 5) then
        if (((status1==0.and.itag==1).and.(mod(it,20)==1)).or.((status1==1.and.itag==1).and.(mod(it,10)==0))) then 
            open(9002,file='p1output.txt',form='formatted',status='unknown',position='append')
                do i = 1,nftnd(1)
                    if (abs(x(1,nsmp(1,i,ift) ))<=38.0d3.and.abs(x(3,nsmp(1,i,ift) )-(-10.0d3))<0.01d0) then
                        write(9002,'(1x,5e32.21e4)') time,fric(71,i,ift),fric(72,i,ift),fric(28,i,ift),fric(29,i,ift)
                    endif
                enddo
            open(9003,file='p2output.txt',form='formatted',status='unknown',position='append')
                do i = 1,nftnd(1)
                    if (abs(x(1,nsmp(1,i,ift) ))<=0.01d0.and.abs(x(3,nsmp(1,i,ift) ))<=40.0d3) then
                        write(9003,'(1x,5e32.21e4)') time,fric(71,i,ift),fric(72,i,ift),fric(28,i,ift),fric(29,i,ift)
                    endif
                enddo                
        endif
    endif

    if (bp == 7) then
        if (((status1==0.and.itag==1).and.(mod(it,20)==1)).or.((status1==1.and.itag==1).and.(mod(it,10)==0))) then 
            open(9002,file='p1output.txt',form='formatted',status='unknown',position='append')
                do i = 1,nftnd(1)
                    if (abs(x(3,nsmp(1,i,ift) )-0.0d3)<0.01d0) then
                        write(9002,'(1x,5e32.21e4)') time,fric(71,i,ift),fric(72,i,ift),fric(28,i,ift),fric(29,i,ift)
                    endif
                enddo
            open(9003,file='p2output.txt',form='formatted',status='unknown',position='append')
                do i = 1,nftnd(1)
                    if (abs(x(1,nsmp(1,i,ift) ))<=0.01d0) then
                        write(9003,'(1x,5e32.21e4)') time,fric(71,i,ift),fric(72,i,ift),fric(28,i,ift),fric(29,i,ift)
                    endif
                enddo                
        endif
    endif
end subroutine output_prof

subroutine output_ruptarea_trac_slip

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j 
    
    if(nftnd(1) > 0) then
        open(unit=1114,file='cplot_ruptarea_trac_slip.txt',status='unknown')    !rupture time
        ! 1,    2,     3,     4
        !fric(81,i),fric(82,i),fric(83,i),fric(84,i)   
        !Ruptured area, total slip, tract at the beginning, tract at the end.
        write(1114,'(1x,4e32.21e4)') (fric(81,i,1),fric(82,i,1),fric(83,i,1),fric(84,i,1), i=1,nftnd(1))    
        close(1114)
    endif
    
end subroutine output_ruptarea_trac_slip

! output_bp8_profiles writes the ten section-4.3 files of SEAS BP8: slip,
! shear stress (both components) and pore pressure along the two cross-section
! lines through the injection point.
!
! Layout, per the benchmark description:
!   row 1        :  0   0   x2(1) x2(2) ... x2(n)
!   rows 2..Nt+1 :  t   log10(max slip rate)   value(1) ... value(n)
subroutine output_bp8_profiles

    use globalvar
    implicit none

    if (bp /= 8) return
    if (nProfRec == 0) return

    call write_one_profile('slip_2_strike.dat',         'Slip_2 (m)',           'slip_2',         'x2', 1, nProfStrike, profS2s)
    call write_one_profile('slip_3_strike.dat',         'Slip_3 (m)',           'slip_3',         'x2', 1, nProfStrike, profS3s)
    call write_one_profile('shear_stress_2_strike.dat', 'Shear_stress_2 (MPa)', 'shear_stress_2', 'x2', 1, nProfStrike, profT2s)
    call write_one_profile('shear_stress_3_strike.dat', 'Shear_stress_3 (MPa)', 'shear_stress_3', 'x2', 1, nProfStrike, profT3s)
    call write_one_profile('pore_pressure_strike.dat',  'Pore_pressure (MPa)',  'pore_pressure',  'x2', 1, nProfStrike, profPs)

    call write_one_profile('slip_2_depth.dat',          'Slip_2 (m)',           'slip_2',         'x3', 2, nProfDepth, profS2d)
    call write_one_profile('slip_3_depth.dat',          'Slip_3 (m)',           'slip_3',         'x3', 2, nProfDepth, profS3d)
    call write_one_profile('shear_stress_2_depth.dat',  'Shear_stress_2 (MPa)', 'shear_stress_2', 'x3', 2, nProfDepth, profT2d)
    call write_one_profile('shear_stress_3_depth.dat',  'Shear_stress_3 (MPa)', 'shear_stress_3', 'x3', 2, nProfDepth, profT3d)
    call write_one_profile('pore_pressure_depth.dat',   'Pore_pressure (MPa)',  'pore_pressure',  'x3', 2, nProfDepth, profPd)

end subroutine output_bp8_profiles

subroutine write_one_profile(fname, colDesc, fieldName, axisName, iline, n, dat)

    use globalvar
    implicit none

    character (len = *), intent(in) :: fname, colDesc, fieldName, axisName
    integer (kind = 4),  intent(in) :: iline           ! 1 = strike line, 2 = dip line
    integer (kind = 4),  intent(in) :: n               ! nodes on that line
    real (kind = dp),    intent(in) :: dat(nProfMax, n)

    integer (kind = 4) :: k, i
    real (kind = dp), allocatable :: coord(:)

    if (iline == 1) then
        allocate(coord(n)); coord = xProfStrike
    else
        allocate(coord(n)); coord = xProfDepth
    endif

    open(9101, file = trim(fname), form = 'formatted', status = 'unknown')
        write(9101,'(A)') '# problem=SEAS Benchmark BP8-QD-'//trim(adjustl(bp8tag))
        write(9101,'(A)') '# author=D.Liu'
        write(9101,'(A)') '# date='//trim(adjustl(runDate))
        write(9101,'(A)') '# code=EQquasi'
        write(9101,'(A)') '# code_version='//EQQUASI_VERSION
        write(9101,'(A,E15.7,A)') '# element_size=', dx, ' m'
        write(9101,'(A)') '# Row #1 = '//trim(axisName)//' (m) with two zeros first'
        write(9101,'(A)') '# Column #1 = Time (s)'
        write(9101,'(A)') '# Column #2 = Max slip rate (log10 m/s)'
        write(9101,'(A,i4,A)') '# Columns #3-', n+2, ' = '//trim(colDesc)
        write(9101,'(A,2(E13.5,A))') '# Computational domain: ', xminc, ' < x2 < ', xmaxc, ' m'
        write(9101,'(A)') '# The line below lists the names of the data fields'
        write(9101,'(A)') trim(axisName)
        write(9101,'(A)') 't'
        write(9101,'(A)') 'max_slip_rate'
        write(9101,'(A)') trim(fieldName)
        write(9101,'(A)') '# Here are the data'

        ! Row 1: two zeros, then the node coordinates.
        write(9101,'(2E22.14,1000E15.7)') 0.0d0, 0.0d0, (coord(i), i = 1, n)
        ! Rows 2..: time, log10 max slip rate, then the field at each node.
        do k = 1, nProfRec
            write(9101,'(E22.14,1000E15.7)') profTime(k), &
                dlog10(max(profVmax(k), 1.0d-30)), (dat(k,i), i = 1, n)
        enddo
    close(9101)

    deallocate(coord)

end subroutine write_one_profile

! output_run_metadata writes runInfo.json next to the results: what was run,
! on what hardware, with how many ranks, and how long it took.
!
! The point is postmortem traceability. A timing quoted without its core count
! and machine load is not reproducible, and a result found on disk months later
! is worth little if nobody can say what produced it.
subroutine output_run_metadata(solverTime, factorTime)

    use globalvar
    implicit none

    real (kind = dp), intent(in) :: solverTime, factorTime
    character (len = 256) :: hostName, cpuModel, cpuLine
    character (len = 64)  :: timeStamp, ompStr
    integer (kind = 4)    :: dtv(8)
    integer (kind = 4)    :: nThreads
    integer (kind = 4)    :: ios
    integer (kind = 4)    :: nCores

    if (me /= 0) return

    call date_and_time(values = dtv)
    write(timeStamp,'(i4.4,A,i2.2,A,i2.2,A,i2.2,A,i2.2,A,i2.2)') &
        dtv(1),'-',dtv(2),'-',dtv(3),'T',dtv(5),':',dtv(6),':',dtv(7)

    hostName = ' '
    cpuModel = ' '
    ! HOSTNAME is not exported to non-interactive shells, so fall back to procfs.
    call get_environment_variable('HOSTNAME', hostName)
    if (len_trim(hostName) == 0) then
        open(unit = 9300, file = '/proc/sys/kernel/hostname', status = 'old', &
             action = 'read', iostat = ios)
        if (ios == 0) then
            read(9300,'(A)',iostat = ios) hostName
            close(9300)
        endif
        if (ios /= 0 .or. len_trim(hostName) == 0) hostName = 'unknown'
    endif

    ! Read the CPU model from /proc/cpuinfo where it exists; leave it unknown
    ! rather than guessing on platforms that have no such file.
    open(unit = 9301, file = '/proc/cpuinfo', status = 'old', action = 'read', iostat = ios)
    if (ios == 0) then
        do
            read(9301,'(A)',iostat = ios) cpuModel
            if (ios /= 0) then
                cpuModel = 'unknown'
                exit
            endif
            if (index(cpuModel, 'model name') > 0) then
                cpuModel = adjustl(cpuModel(index(cpuModel,':')+1:))
                exit
            endif
        enddo
        close(9301)
    else
        cpuModel = 'unknown'
    endif

    ! Report the requested OpenMP thread count rather than linking against the
    ! runtime, so this stays correct whether or not the build enabled OpenMP.
    nCores = 0
    open(unit = 9303, file = '/proc/cpuinfo', status = 'old', action = 'read', iostat = ios)
    if (ios == 0) then
        do
            read(9303,'(A)',iostat = ios) cpuLine
            if (ios /= 0) exit
            if (index(cpuLine, 'processor') == 1) nCores = nCores + 1
        enddo
        close(9303)
    endif

    ! 0 means OMP_NUM_THREADS was not set, i.e. each rank may spawn as many
    ! threads as it likes. That is not the same as 1, and on a shared node it is
    ! how a run ends up oversubscribed -- worth recording as its own value
    ! rather than silently reporting 1.
    nThreads = 0
    ompStr = ' '
    call get_environment_variable('OMP_NUM_THREADS', ompStr)
    if (len_trim(ompStr) > 0) then
        read(ompStr,*,iostat = ios) nThreads
        if (ios /= 0) nThreads = 0
    endif

    open(9302, file = 'runInfo.json', form = 'formatted', status = 'unknown')
        write(9302,'(A)') '{'
        write(9302,'(A)')      '  "code": "EQquasi",'
        write(9302,'(A)')      '  "version": "'//EQQUASI_VERSION//'",'
        write(9302,'(A,I0,A)') '  "benchmark_id": ', bp, ','
        write(9302,'(A)')      '  "run_timestamp": "'//trim(timeStamp)//'",'
        write(9302,'(A)')      '  "host": "'//trim(hostName)//'",'
        write(9302,'(A)')      '  "cpu_model": "'//trim(adjustl(cpuModel))//'",'
        write(9302,'(A,I0,A)') '  "host_logical_cores": ', nCores, ','
        write(9302,'(A,I0,A)') '  "mpi_ranks": ', nprocs, ','
        write(9302,'(A,I0,A)') '  "omp_threads_per_rank": ', nThreads, ','
        write(9302,'(A,I0,A)') '  "num_nodes": ', numnp, ','
        write(9302,'(A,I0,A)') '  "num_elements": ', numel, ','
        write(9302,'(A,I0,A)') '  "num_fault_nodes": ', nftnd(1), ','
        write(9302,'(A,I0,A)') '  "num_equations": ', neq, ','
        write(9302,'(A,E15.7,A)') '  "element_size_m": ', dx, ','
        write(9302,'(A,I0,A)') '  "steps_completed": ', it-1, ','
        write(9302,'(A,E15.7,A)') '  "simulated_time_s": ', time, ','
        write(9302,'(A,E15.7,A)') '  "time_loop_seconds": ', solverTime, ','
        write(9302,'(A,E15.7,A)') '  "factorization_seconds": ', factorTime, ','
        write(9302,'(A,E15.7,A)') '  "seconds_per_step": ', solverTime/max(1,it-1), ','
        write(9302,'(A,E15.7,A)') '  "max_slip_rate_final_m_s": ', maxSlipRate, ','
        write(9302,'(A,I0,A)') '  "solver": ', sol_op, ','
        write(9302,'(A,I0,A)') '  "friclaw": ', friclaw, ','
        write(9302,'(A,I0)')   '  "fluid_source_model": ', fluid_src
        write(9302,'(A)') '}'
    close(9302)

    write(*,*) '=     Run metadata written to runInfo.json                          ='

end subroutine output_run_metadata
