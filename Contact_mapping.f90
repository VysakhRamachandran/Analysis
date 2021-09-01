       program contact_map

       implicit none

       real(kind=8),allocatable :: rMG(:,:,:),rRNA(:,:,:)
       real(kind=8):: boxl1(3),boxl2(3),dr(3),cutoff1,cutoff2,dist,summ=0.0,actual_contact_no,theta,z_projection,xy_projection
       
       
       integer :: natom_MG,natom_RNA,iatom_MG,iatom_RNA,count1,av_contact
       integer :: iframe,nframe,i
      
     
       character(len=4) :: junk2
       character(len=6) :: junk1

       open(unit=2,file='Outer_Sphere_MG97_0_255.gro', action='read')
       open(unit=3,file='Outer_Sphere_All_RNA_Atoms0_255.gro', action='read')
       !open(unit=5,file='Density_Map_MG96_0_134.dat', action='write', status='old')
       !open(unit=10,file='Density_Map_MG96_0_134_atom1.dat', action='write', status='old')
       !open(unit=11,file='Density_Map_MG96_0_134_atom2.dat', action='write', status='old')
       !open(unit=20,file='VMD_DOS_MG97_All_RNA_Atoms_FM_0_255nm_rdfShell0_37nm_to_0_55nm.dat', action='write', status='new')
       !open(unit=25,file='2D_DOS_MG97_All_RNA_Atoms_FM_0_255nm_rdfShell0_37nm_to_0_55nm.dat', action='write', status='new')
       open(unit=50, file='MG97_All_RNA_Atoms_ND_FM_0_255nm_rdfShell0_18nm_to_0_22nm.dat', action='write', status='new')   
    
       ! System Information
       ! Number of MG Atoms
       natom_MG= 1
       natom_RNA= 614
       write(6,*)'MG Atom',natom_MG
       write(6,*)'RNA Atom',natom_RNA


       !Total Number of Snapshot
       nframe= 50
       write(6,*)'Total Frame',nframe
       cutoff1=1.8d0
       cutoff2=2.2d0
       !cutoff3=3.0d0
       !cutoff4=3.6d0
       !cutoff5=3.7d0
       !cutoff6=5.5d0
       !cutoff7=5.6d0
       !cutoff8=8.0d0

       ! Allocation
        allocate(rMG(natom_MG,nframe,3),rRNA(natom_RNA,nframe,3))
        rMG=0.0d0;rRNA=0.0d0

     

       
        do iframe = 1 , nframe
       ! Reading MG Coordinates
        read(2,*)
        read(2,*)
        do iatom_MG = 1 , natom_MG
        read(2,*) junk1, junk2, i, rMG(iatom_MG,iframe,:)
        rMG(iatom_MG,iframe,:)=rMG(iatom_MG,iframe,:)*10.0d0
        end do
        read(2,*)boxl1(:)
        
        ! Converting into \AA
        boxl1=boxl1*10.0d0

        ! Reading RNA atom's Coordinates
        read(3,*)
        read(3,*)
        do iatom_RNA = 1 , natom_RNA
        read(3,*) junk1, junk2, i, rRNA(iatom_RNA,iframe,:)
        rRNA(iatom_RNA,iframe,:)=rRNA(iatom_RNA,iframe,:)*10.0d0
        end do
        read(3,*)boxl2(:)
        end do
        ! Converting into \AA
        boxl2=boxl2*10.0d0

          
          do iframe =1, nframe
             count1=0
             do iatom_MG = 1 ,natom_MG
                do iatom_RNA = 1 , natom_RNA
                   dr(:)= rMG(iatom_MG,iframe,:)-rRNA(iatom_RNA,iframe,:)
                   dr(:)= dr(:)-boxl1(:)*anint(dr(:)/boxl1(:))
                   dist= sqrt(dot_product(dr,dr))
                   if(dist .ge. cutoff1 .and. dist .le. cutoff2)then
                     !print*, dr(:)
                     theta=acos(dr(3)/dist)
                     z_projection=dist*cos(theta)      
                     xy_projection=dist*sin(theta)            
                   if (iframe==nframe) then
                     !write(20,*) iatom_RNA
                     !write(25,*) xy_projection,z_projection
                   end if
                     !write(5,*) rRNA(iatom_RNA,iframe,:), dist
                     !print*, rRNA(iatom_RNA,iframe,:),iatom_RNA,dr(:)
                     count1=count1+1
                   !if (count1==1) then 
                     !write(10,*) dr(:) !,dist,rRNA(iatom_RNA,iframe,:),iatom_RNA
                   !elseif (count1==2) then
                     !write(11,*) dr(:) !,dist,rRNA(iatom_RNA,iframe,:),iatom_RNA
                   !end if
                   endif
                
                enddo
             end do
         print*,count1
         write(50,*) iframe,count1
         summ=summ+count1   
         !print*,summ      
         end do
         actual_contact_no=(summ/nframe)
         av_contact=NINT((summ/nframe))  
         print*,"So the actual no. of contacts after doing the average are",actual_contact_no
         print*,"Avrage no. of contacts taking nearest integer of 'actual_contact_no' are",av_contact
       end program contact_map
