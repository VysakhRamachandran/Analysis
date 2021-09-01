       program Ct_St_TCF

       implicit none

       real(kind=8),allocatable :: rprot(:,:),rw(:,:,:),dist(:)
       real(kind=8),allocatable :: correlation(:)
       real(kind=8):: boxl(3),boxl1(3),boxl2(3),pi,dr(3),cutoff,av
       real(kind=8):: inter_oxy_dist,dist11,dist12,dist21,dist22
       real(kind=8):: temp_array(4),angle_conv,cond_3,ohdistsq
       real(kind=8):: u1(3),u2(3),min_angle,oodistsq,min_dist,sum2
       real(kind=8):: corr,sum1,dt=2.0d0,p,q,cutoff1,cont,summ=0.0,Number_bulk_water

       integer(kind=1),allocatable :: location(:,:)
       integer,allocatable :: ht_val(:,:)
       integer,allocatable :: info(:,:,:)!,hb_info(:,:)
       integer :: i1, i2, gw, rgw,jmol_water,corr_frame,ip,j,m,count2
       integer :: natom_prot,iatom_prot,natom_water,iatom_water,k
       integer :: nmol_water,imol_water,iframe,nframe,i,ii,pair,jframe
       integer :: count1,count4
      
       character(len=4),allocatable :: prot_atom(:),acceptor(:)
       character(len=4) :: junk2
       character(len=6) :: junk1
       character(len=100) :: prot_gro,sol_gro
      
       open(unit=2,file='All_RNA_atoms_MG97_conf0_261.gro', action='read') 
       open(unit=3,file='GammaCal_water_MG97_conf0_261.gro', action='read')
       open(unit=10, file='Number_Bulk_water_MG97_FM.dat', action='write', status='new')
       ! System Information
       ! Number of Protein Atoms (Lysozyme)
       natom_prot = 3047
       write(69,*)'Protein Atom',natom_prot
       ! Number of Water molecules
       nmol_water=31711;natom_water=3
       write(69,*)'Water Molecule',nmol_water
       ! Total Number of Snapshot
       nframe= 1000
       write(69,*)'Total Frame',nframe
       ! Constants
       cutoff=20.0d0

       ! Allocation
       allocate(rw(nmol_water,natom_water,3), rprot(natom_prot,3))
       allocate(prot_atom(natom_prot),location(nframe,nmol_water))
       allocate(dist(natom_prot))
       rw=0.0d0; rprot=0.0d0; prot_atom='zz';
       
       ! Reading the Protein atom name from 1hew.gro
       !read(3,*)
       !read(3,*)
       !do iatom_prot = 1 , natom_prot
        !read(3,*) junk1, prot_atom(iatom_prot)
       !end do
       ! ==========================
frm:   do iframe = 1 , nframe
          count1=0
        ! Reading Protein Coordinates
        read(2,*)
        read(2,*)
        do iatom_prot = 1 , natom_prot
         read(2,*) junk1,junk2,i,rprot(iatom_prot,:)
        enddo
        read(2,*)boxl1(:)
        ! Converting into \AA
        rprot=rprot*10.0d0;boxl1=boxl1*10.0d0
        ! Reading Water coordinates
        read(3,*)
        read(3,*)
        do imol_water = 1 , nmol_water
         do iatom_water = 1 , natom_water
          if(imol_water <= 3333) then
          read(3,*) junk1, junk2, i, rw(imol_water,iatom_water,:)
          else
          read(3,*) junk1, junk2,  rw(imol_water,iatom_water,:)
          endif
         end do
        end do
        read(3,*)boxl2(:)
        ! Converting into \AA
        rw=rw*10.0d0;boxl2=boxl2*10.d0
        if(mod(iframe,500)==0)then 
         write(69,*)'Frame No.',iframe
         write(69,*)boxl1
         write(69,*)boxl2
        end if
        ! Searching the water 
         do imol_water = 1 ,nmol_water
          do iatom_water = 1 ,natom_water
             dist=999.99d0
          do iatom_prot = 1 , natom_prot
           dr(:)=rw(imol_water,1,:)-rprot(iatom_prot,:)
           dr(:)= dr(:)-boxl1(:)*anint(dr(:)/boxl1(:))
           dist(iatom_prot)= sqrt(dot_product(dr,dr))
          enddo
          enddo
         if(minval(dist)>= cutoff)then
           count1=count1+1
         end if
         enddo
         write(10,*) iframe,count1
         summ=summ+count1   
         end do frm
         Number_bulk_water=(summ/nframe) 
         write(10,*) Number_bulk_water
         print*,"So the actual no. of bulk water after doing the average is",Number_bulk_water
       end program Ct_St_TCF  
