subroutine output_onfault_st

    use globalvar
    implicit none
    
    integer (kind = 4) :: i, j 

    if(n4onf > 0) then
        do i=1,n4onf
            j=anonfs(3,i)
            if(j==1)  then  !main fault stations
                sttmp = '      '
                dptmp = '      '
                                if (bp == 7) then 
                  write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)) 
                  write(dptmp,'(i4.3)') int(-xonfs(2,anonfs(2,i),j)) 
                                else
                                  write(sttmp,'(i4.3)') int(xonfs(1,anonfs(2,i),j)/1000.d0)
                                  write(dptmp,'(i4.3)') int((-xonfs(2,anonfs(2,i),j))/1000.d0)
                                endif
                open(51,file='fltst_strk'//trim(adjustl(sttmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

                ! sttmp = '      '
                ! dptmp = '      '
                ! write(sttmp,'(f5.1)') xonfs(1,anonfs(2,i),j)/1000.d0 
                ! write(dptmp,'(f5.1)') abs(xonfs(2,anonfs(2,i),j)/1000.d0) 
                ! loca = '# location = on fault, '//trim(adjustl(sttmp))//' km along strike, '//trim(adjustl(dptmp))//' km down-dip'        
            endif
                        !write(51,*) '# This is the file header'
            !write(51,*) '# problem=SEAS Benchmark No.4'
            ! write(51,*) '#Author=D.Liu'
                !write(51,*) '# code = EQSimu'
            !write(51,*) '# version=1.0'
            !write(51,*) '# modeler = A.Modeler'
                        !write(51,*) '#element_size = 200 m x 200 m on fault'
            ! write(51,*) '# Time series in 8 columns in format e15.7'
            ! write(51,*) '# Column #1 = Time (s)'
            ! write(51,*) '# Column #2 = horizontal slip (m)'
            ! write(51,*) '# Column #3 = horizontal slip rate (m/s)'
            ! write(51,*) '# Column #4 = horizontal shear stress (MPa)'
            ! write(51,*) '# Column #3 = down-dip slip (m)'
            ! write(51,*) '# Column #6 = down-dip slip rate (m/s)'
            ! write(51,*) '# Column #7 = down-dip shear stress (MPa)'
            ! !write(51,*) '# Column #8 = normal stress (MPa)'
            ! write(51,*) '# The line below lists the names of the data fields:'
            ! write(51,*) 't h-slip h-slip-rate h-shear-stress v-slip v-slip-rate v-shear-stress'
            ! write(51,*) '#'        
            do j=1,it-1
            write(51,'( 8e32.21e4)') fltsta(1,j,i),fltsta(5,j,i),-fltsta(6,j,i), fltsta(2,j,i), &
                                fltsta(3,j,i),fltsta(8,j,i)/1.d6,-fltsta(9,j,i)/1.d6,dlog10(fltsta(4,j,i))
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
                        if (bp == 7) then
              write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i))) 
              write(sttmp,'(i4.3)')   int(x4nds(1,an4nds(1,i))) 
              write(dptmp,'(i4.3)')   int(abs(x4nds(3,an4nds(1,i)))) 
                        else
                          write(bodytmp,'(i4.3)') int(x4nds(2,an4nds(1,i))/1000.d0)
                          write(sttmp,'(i4.3)')   int(x4nds(1,an4nds(1,i))/1000.d0)
                          write(dptmp,'(i4.3)')   int(x4nds(3,an4nds(1,i))/1000.d0)
                        endif
            
            open(51,file='srfst_strk'//trim(adjustl(sttmp))//'st'//trim(adjustl(bodytmp))//'dp'//trim(adjustl(dptmp))//'.txt',status='unknown')

            ! bodytmp = '      '
            ! sttmp = '      '
            ! dptmp = '      '
            ! write(bodytmp,'(f5.1)') x4nds(2,an4nds(1,i))/1000. 
            ! write(sttmp,'(f5.1)') x4nds(1,an4nds(1,i))/1000. 
            ! write(dptmp,'(f5.1)') abs(x4nds(3,an4nds(1,i)))/1000. 
            ! loca = '# location = '//trim(adjustl(bodytmp))//' km off fault, '//trim(adjustl(sttmp))//' km along strike'//trim(adjustl(dptmp))//' km depth'
            ! write(51,*) '# ',projectname
            ! write(51,*) '# Author=',author
            ! call date_and_time(values=time_array)
            ! write(51,'( a10,i2,a1,i2,a1,i4,a1,i2,a1,i2,a1,i2)') ' # date = ',time_array(2), &
                    ! '/',time_array(3),'/',time_array(1),' ',time_array(5),':',time_array(6), &
                    ! ':',time_array(7)
            ! write(51,*) '# code = EQdyna3D'
            ! write(51,*) '# code_version = 5.1.0'
            ! write(51,*) '# element_size =',dx
            ! write(51,'( a14,f8.4,a3)') '# time_step=', dt, ' s'
            ! write(51,'( a19,i6)') '# num_time_steps=',locplt
            ! !write(51,*) loca
            ! !write(51,*) '# Time series in 11 columns in format E15.7 for data'
            ! write(51,*) '# Column #1 = Time (s)'
            ! write(51,*) '# Column #2 = horizontal displacement (m)'
            ! write(51,*) '# Column #3 = horizontal velocity (m/s)'
            ! write(51,*) '# Column #4 = vertical displacement (m)'
            ! write(51,*) '# Column #5 = vertical velocity (m/s)'
            ! write(51,*) '# Column #6 = normal displacement (m)'
            ! write(51,*) '# Column #7 = normal velocity (m/s)'
            ! write(51,*) '#'
            ! write(51,*) '# The line below lists the names of the data fields:'
            ! write(51,*) 't h-disp h-vel v-disp v-vel n-disp n-vel'
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
    
    integer (kind = 4) :: i, j 
    
    if(nftnd(1) > 0) then
        open(unit=1111,file='cplot_EQquasi.txt',status='unknown')    !rupture time
        ! write(1111,*) '#This is the file header:'
        ! write(1111,*) '# problem=TPV8'
        ! write(1111,*) '# author=D.Liu'
        ! write(1111,*) '# date=2016/03/10'
        ! write(1111,*) '# code=EQquasi3d'
        ! write(1111,*) '# code_version=1.0'
        ! write(1111,*) '# element_size=100 m'
        ! write(1111,*) '# Column #1 = horizontal coordinate, distance along strike (m)'
        ! write(1111,*) '# Column #2 = vertical coordinate, distance down-dip (m)'
        ! write(1111,*) '# Column #3 = rupture time (s)'
        ! write(1111,*) '# The line below lists the names of the data fields.'
        ! write(1111,*) '# It indicates that the first column contains the horizontal'
        ! write(1111,*) '# coordinate (j), the second column contains the vertical'
        ! write(1111,*) '# coordinate (k), and the third column contains the time (t).'
        ! write(1111,*) 'j k t'
        ! write(1111,*) '# Here is the rupture history'    
        !write(1111,'(1x,7e32.21e4)') ((x(j,nsmp(1,i,1)),j=1,3),fnft(i,1),&
        !    (fltslp(j,i,1),j=1,3),i=1,nftnd(1))
        !!-FORMAT
        !!-1,    2,     3, fric(26,i),4,fric(20,i),5,fric(28,i),6,fric(29,i),7,fric(30,i)   
        !!-XCOOR,+ZCOOR, V_trial,     psi,        ,tstk,       ,tdip,       ,tnrm0
        write(1111,'(1x,16e32.21e4)') (x(1,nsmp(1,i,1)),-x(3,nsmp(1,i,1)),fric(26,i,1),fric(20,i,1),fric(28,i,1), &
                fric(29,i,1),fric(30,i,1),fric(31,i,1),fric(32,i,1),fric(33,i,1),fric(34,i,1),fric(35,i,1),fric(36,i,1), &
                                fric(44,i,1), fric(45,i,1), fnft(i,1), i=1,nftnd(1))    
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
    
    open(1113,file='global.dat',form='formatted',status='unknown')
            write(1113,'(7e32.21e4)') (globaldat(1,i), globaldat(2,i), globaldat(3,i), globaldat(4,i), globaldat(5,i), &
                                           globaldat(6,i), globaldat(7,i), i=1,it-1)
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
