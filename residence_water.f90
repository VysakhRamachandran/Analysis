      program Ct_St_TCF

       implicit none

       real(kind=8),allocatable :: rprot(:,:),rw(:,:,:),dist(:),dist1(:)
       real(kind=8),allocatable :: counter(:),correlation(:)
       real(kind=8):: boxl(3),boxl1(3),boxl2(3),pi,dr(3),cutoff
       real(kind=8):: inter_oxy_dist,dist11,dist12,dist21,dist22
       real(kind=8):: temp_array(4),angle_conv,cond_3,ohdistsq
       real(kind=8):: u1(3),u2(3),min_angle,oodistsq,min_dist,sum2
       real(kind=8):: corr,sum1,dt=20.0d0,p,q

       integer,allocatable :: gwi(:),rgwi(:),stay(:)
       integer(kind=1),allocatable :: location(:,:)
       integer,allocatable :: ht_val(:,:)
       integer,allocatable :: info(:,:,:)!,hb_info(:,:)
       integer :: i1, i2, gw, rgw,jmol_water,corr_frame,ip,j,m,count2,count1
       integer :: natom_prot,iatom_prot,natom_water,iatom_water,k
       integer :: nmol_water,imol_water,iframe,nframe,i,ii,pair,jframe
       integer :: hbpair,kk,ih,flag
      
       character(len=4),allocatable :: prot_atom(:),acceptor(:)
       character(len=4) :: junk2
       character(len=6) :: junk1
       character(len=100) :: prot_gro,sol_gro
      
       open(unit=2,file='P.gro', action='read') 
       open(unit=3,file='K.gro', action='read')
       
       ! System Information
       ! Number of Protein Atoms (Lysozyme)
       natom_prot = 93
       write(67,*)'Protein Atom',natom_prot
       ! Number of Water molecules
       nmol_water=138;natom_water=1
       write(67,*)'Water Molecule',nmol_water
       ! Total Number of Snapshot
       nframe= 9000
       write(67,*)'Total Frame',nframe
       ! Constants
       pi=4.0d0*atan(1.0d0)
       angle_conv=180.0d0/pi
       cutoff=4.8d0
       oodistsq=3.5d0*3.5d0
       ohdistsq=2.45d0*2.45d0
       cond_3=30.0d0

       ! Allocation
       allocate(rw(nmol_water,natom_water,3), rprot(natom_prot,3))
       allocate(prot_atom(natom_prot),location(nframe,nmol_water))
       allocate(dist(natom_prot))
       allocate(stay(nframe))
       allocate(gwi(nmol_water),rgwi(nmol_water))
       rw=0.0d0; rprot=0.0d0; prot_atom='zz';location=0
       gwi=0;rgwi=0
       ! Reading the Protein atom name from 1hew.gro
       !read(3,*)
       !read(3,*)
       !do iatom_prot = 1 , natom_prot
        !read(3,*) junk1, prot_atom(iatom_prot)
       !end do
       ! ==========================
frm:   do iframe = 1 , nframe
        ! Reading Protein Coordinates
        read(2,*)
        read(2,*)
        do iatom_prot = 1 , natom_prot
         read(2,*) junk1, junk2, i, rprot(iatom_prot,:)
        end do
        read(2,*)boxl1(:)
        ! Converting into \AA
        rprot=rprot*10.0d0;boxl1=boxl1*10.0d0
        ! Reading Water coordinates
        read(3,*)
        read(3,*)
        do imol_water = 1 , nmol_water
         do iatom_water = 1 , natom_water
          read(3,*) junk1, junk2, i, rw(imol_water,iatom_water,:)
         end do
        end do
        read(3,*)boxl2(:)
        ! Converting into \AA
        rw=rw*10.0d0;boxl2=boxl2*10.d0
        if(mod(iframe,500)==0)then 
         write(67,*)'Frame No.',iframe
         write(67,*)boxl1
         write(67,*)boxl2
        end if
        ! Searching the water 
        gw = 0; gwi=0
        do imol_water = 1 ,nmol_water
         dist=999.99d0
         do iatom_prot = 1 , natom_prot
          dr(:)=rw(imol_water,1,:)-rprot(iatom_prot,:)
          dr(:)= dr(:)-boxl1(:)*anint(dr(:)/boxl1(:))
          dist(iatom_prot) = sqrt(dot_product(dr,dr))
         end do
         if(minval(dist) < cutoff)then
          gw = gw + 1
          gwi(gw) = imol_water
          location(iframe,imol_water)=1
         end if
        end do
        if(iframe == 1)then
         write(67,*)'At t=0, Water within dist=',gw
         rgw=gw
         rgwi(1:gw)=gwi(1:gw)
        else
         do imol_water = 1 , gw
          if(any(rgwi==gwi(imol_water)))cycle
          rgw = rgw + 1
          rgwi(rgw)=gwi(imol_water) 
         end do
        end if
        !write(*,*)gw
       end do frm
       write(67,*)'Total # Water molecules over the trajectory',rgw
       do i=1,rgw
        write(210,*) i, rgwi(i)
       enddo
       do i=1,nframe
         do j=1,nmol_water
          write(211,*) location(i,j)
         enddo
       enddo
     
         
         do ip=1,rgw
              imol_water=rgwi(ip)
                 stay=0
!                k=1;count2=0 You have to start with k=0, because in the next loop it should increase one-by-one. see next
                  k=0; count2 =0
               do iframe=1,nframe
                  k = k+1 !you have to add this line to increase the fr
!                 if((iframe+jframe)>nframe+1) exit  ! this line has no meaning as jframe information is not there, as it is starting next loop, see the next line
                  if((k+1)>nframe) exit 
!               if(location(iframe,imol_water)==0)cycle ! this line can't work, becoz you have no info about jframe as it is starting from next loop, it should be modified as below
               if(location(k,imol_water)==0)cycle
               count1=count1+1;count2=count2+1
              do jframe = k+1,nframe
               if(location(jframe,imol_water)==1)then 
!               count1=count1+1
               stay(count2)=stay(count2) + 1
!              stay(count2)=count1, you can not have this line here, bocoz this will come after final counting of a stay, see below
              else
              k = jframe !you have to add this line, to continue the next loop, just check
!              k=jframe+1
              exit 
              end if
              k = jframe
              end do
              end do 
              write(212,*) imol_water, maxval(stay)  
              write(213,*) maxval(stay)
               p = maxval(stay)
               if (p> 200) then    
               write(217,*) imol_water, maxval(stay)
               write(218,*) maxval(stay)
               endif 
               if ( p <= 200) then
               write(219,*) imol_water, maxval(stay)
               write(220,*) maxval(stay)
               endif 
                100 format ( (i5,a10,i5,3f8.3))     
              end do
                  
              
            
            
                
             end program Ct_St_TCF  
