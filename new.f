c-----------------------------------------------------------------------
C
C  USER SPECIFIED ROUTINES:
C
C     - boundary conditions
C     - initial conditions
C     - variable properties
C     - local acceleration for fluid (a)
C     - forcing function for passive scalar (q)
C     - general purpose routine for checking errors etc.
C
c-----------------------------------------------------------------------
      subroutine uservp (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      udiff =0.
      utrans=0.
      return
      end
c-----------------------------------------------------------------------
      subroutine userf  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)


c     Note: this is an acceleration term, NOT a force!
c     Thus, ffx will subsequently be multiplied by rho(x,t).


      real retau
c
      visc=param(2)
      g=param(5)
      retau=310
      ffx = 0.0
      ffy = 0.0
      ffz = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userq  (ix,iy,iz,eg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      qvol   = 0.0
      source = 0.0

      return
      end
c-----------------------------------------------------------------------
      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
c      common /rparts/ pts(ldim,lpart),vel(ldim,2:3,lpart)
c      common /iparts/ npart,partid(lpart)
      real   xerange(2,3,lelt)
      common /elementrange/ xerange

      if(istep.eq.0)then
        ntot = lx1*ly1*lz1*nelt
        nxyz = lx1*ly1*lz1
        do ie = 1,nelt
           xerange(1,1,ie) = vlmin(xm1(1,1,1,ie),nxyz)
           xerange(2,1,ie) = vlmax(xm1(1,1,1,ie),nxyz)
           xerange(1,2,ie) = vlmin(ym1(1,1,1,ie),nxyz)
           xerange(2,2,ie) = vlmax(ym1(1,1,1,ie),nxyz)
           xerange(1,3,ie) = vlmin(zm1(1,1,1,ie),nxyz)
           xerange(2,3,ie) = vlmax(zm1(1,1,1,ie),nxyz)
        enddo  
      endif

      
      
       !call concentration 
       !call exitt

      ifxyo = .true.
      if (istep.gt.iostep) ifxyo = .false.
      if (istep.eq.0) call outpost(vx,vy,vz,pr,t,'   ') ! Part. coordination
      
      call my_particle_generator                ! Particle injection
      


      return
      end
c-----------------------------------------------------------------------
      subroutine userbc (ix,iy,iz,iside,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
            
      
      one = 1.0
      pi  = 4.*atan(one)
      uy  = (0.25)*(pi*pi)*sin(pi*(0.5-x))*sin(pi*(0.5-z))
      ux  = 0.0
      uz  = 0.0

      temp=0
      

      return
      end
c-----------------------------------------------------------------------
      subroutine useric (ix,iy,iz,ieg)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      one = 1.0
      pi  = 4.*atan(one)
      uy  = (0.25)*(pi*pi)*sin(pi*(0.5-x))*sin(pi*(0.5-z))
      ux=0.0
      
      uz=0.0
      temp=0.0
      
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat
      include 'SIZE'
      include 'TOTAL'
      
      param(66) = 0.
      param(67) = 0.
c
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'
c
      return
      end
c-----------------------------------------------------------------------
      subroutine usrdat3  !  Modify geometry 

      include 'SIZE'
      include 'TOTAL'


      return
      end
c----------------------------------------------------------------------
      subroutine set_part_pointers
      include 'SIZE'
c     Minimal value of ni = 5
c     Minimal value of nr = 14*ndim + 1
      common  /iparti/ n,nr,ni
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      jrc = 1 ! Pointer to findpts return code
      jpt = 2 ! Pointer to findpts return processor id
      je0 = 3 ! Pointer to findpts return element id
      jps = 4 ! Pointer to proc id for data swap
      jai = 5 ! Pointer to auxiliary integers
      nai = ni - (jai-1)  ! Number of auxiliary integers
      if (nai.le.0) call exitti('Error in nai:$',ni)

      jr  = 1         ! Pointer to findpts return rst variables
      jd  = jr + ndim ! Pointer to findpts return distance
      jx  = jd + 1    ! Pointer to findpts input x value
      jy  = jx + 1    ! Pointer to findpts input y value
      jz  = jy + 1    ! Pointer to findpts input z value

      jx1 = jx + ndim ! Pointer to xyz at t^{n-1}
      jx2 = jx1+ ndim ! Pointer to xyz at t^{n-2}
      jx3 = jx2+ ndim ! Pointer to xyz at t^{n-3}
      
      ja0 = jx3+ ndim ! Pointer to current particle acceleration
      ja1 = ja0+ ndim ! Pointer to particle acceleration at t^{n-1}
      ja2 = ja1+ ndim ! Pointer to particle acceleration at t^{n-2}
      ja3 = ja2+ ndim ! Pointer to particle acceleration at t^{n-3}
      
      jv0 = ja3+ ndim ! Pointer to current particle velocity
      jv1 = jv0+ ndim ! Pointer to particle velocity at t^{n-1}
      jv2 = jv1+ ndim ! Pointer to particle velocity at t^{n-2}
      jv3 = jv2+ ndim ! Pointer to particle velocity at t^{n-3}

      ju0 = jv3+ ndim ! Pointer to current fluid velocity
      ju1 = ju0+ ndim ! Pointer to fluid velocity at t^{n-1}
      ju2 = ju1+ ndim ! Pointer to fluid velocity at t^{n-2}
      ju3 = ju2+ ndim ! Pointer to fluid velocity at t^{n-3}

      jf0 = ju3+ ndim ! Pointer to forcing at current timestep

      jar = jf0+1 ! Pointer to auxiliary reals
      nar = nr - (jar-1)  ! Number of auxiliary reals

      if (nar.le.0) call exitti('Error in nar:$',nr)
      return
      end
c----------------------------------------------------------------------
	  subroutine concentration (rpart,nr,ipart,ni,n)
	  include 'SIZE'
      include 'TOTAL'
      real    rpart(nr,n)
      integer ipart(ni,n)
      real   xerange(2,3,lelt)
      common /elementrange/ xerange
      real concen(lelt),xcon(lx1,ly1,lz1,lelt)
      real ycon(lx1,ly1,lz1,lelt), zcon(lx1,ly1,lz1,lelt)
      common /concon/ concen,xcon,ycon,zcon
      
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      integer ie
      real vp,pi,pd,total_con
      pi    = 4.*atan(1.0)
      pd = .01
      vp=(4/3)*pi*(pd/2)**3
      total_con=(vp*n)/(volvm1-vp*n)
      if(istep.eq.0) then 
      if (nid.eq.0) write(6,*)'total concentration is: ', total_con
      endif
      !if(icalld.eq.0) rzero(concen,nelt)
      write(6,*) 'made it here!!' 
      do ip = 1,n
         xloc = rpart(jx,ip)
         yloc = rpart(jy,ip)
         zloc = rpart(jz,ip)
         do ie=1,nelt
            if (xloc.ge.xerange(1,1,ie).and.xloc.le.xerange(2,1,ie))then
            if (yloc.ge.xerange(1,2,ie).and.yloc.le.xerange(2,2,ie))then
            if (zloc.ge.xerange(1,3,ie).and.zloc.le.xerange(2,3,ie))then
         concen(ie)=concen(ie)+log(vp/(volel(ie)*dt*(1/total_con)))
                write(6,'(A,i5,A,e12.7)') 'concentration in element: '
     &          ,ie,' is ',concen(ie)                 
                goto 123
            endif
            endif
            endif
         enddo
123      continue
      enddo
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt
c      ipstep = iostep/(iostep/100)
c      if(mod(istep,ipstep).eq.0.or.istep.eq.0)  then 
      do ie=1,nelt
      do i=1,nxyz
            xcon(i,1,1,ie)=concen(ie)
            ycon(i,1,1,ie)=concen(ie)
            zcon(i,1,1,ie)=concen(ie)
      enddo
      enddo
         call outpost (xcon,ycon,zcon,pr,t,      'con')  
c      endif       
c      do e=1,nelt
c         vv=vv+volel(e)
c         write(6,'(e12.6,I6,A)') volel(e) ,e, ' volume of element'
c      enddo
c         write(6,'(A,e12.6)') 'total volume is: ',vv  
	  return
	  end
	  
c----------------------------------------------------------------------
      subroutine particle_find(rpart,nr,ipart,ni,n)
      include 'SIZE'
      include 'TOTAL'
      
      real    rpart(nr,n)
      integer ipart(ni,n)
      parameter (lrf=18*ldim+2,lif=6)
c      real               rfpts(lrf,lpart)
c      common /fptspartr/ rfpts
c      integer            ifpts(lif,lpart),fptsmap(lpart)
c      common /fptsparti/ ifpts,fptsmap
      
      parameter(nmax=lpart) 
      common /iv_intp/ ihandle 
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      integer icalld
      save    icalld
      data    icalld /0/
            
      
c      write(6,*) jrc ,'jcr'
c      write(6,*) jx, 'jx'
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt
     
      if (icalld.eq.0) then		! interpolation setup	!? intpts_done(ih_intp_v)?
        icalld = 1
        tolin  = 1.e-12
        if (wdsize.eq.4) tolin = 1.e-6
        call intpts_setup(tolin,ihandle)
      endif
c       write(6,*) lrf,'lrf'
c       write(6,*) lif,'lif'
c      write(6,*) 'made it here 2'
        if(nid.eq.0) write(6,*) 'call findpts'
         call findpts(ihandle,ipart(jrc,1),lif,
     &               ipart(jpt,1),lif,
     &               ipart(je0,1),lif,
     &               rpart(jr, 1),lrf,
     &               rpart(jd, 1),lrf,
     &               rpart(jx, 1),lrf,
     &               rpart(jy, 1),lrf,
     &               rpart(jz, 1),lrf, n)
     
        nfail = 0 
       do in=1,n
           ! check return code
c           if(ipart(jrc,in).eq.1) then
c             if(rpart(jd,in).gt.1e-12) then 
c               nfail = nfail + 1
c               write(6,'(1p3e15.7,I2)') 
c     &     ' WARNING: point on boundary or outside the mesh xy[z]d^2: ',
c     &     (rpart(jx+k,in),k=0,ndim-1),ipart(jrc,in)
c             endif   
           
        enddo
     
         
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt

      if (n.gt.nmax) call exitti ('ABORT: interp_v() n > nmax!$',n)
      
      if (nelgt.ne.nelgv) call exitti
     $   ('ABORT: interp_v() nelgt.ne.nelgv not yet supported!$',nelgv) 
     
     
      return
      end
c-------------------------------------------------------------------------      
      subroutine my_particle_generator ! Particle injection
      include 'SIZE'
      include 'TOTAL'
      include 'mpif.h'
      
      parameter (lr=18*ldim+2,li=5+1)
      common /rparts/ rpart(lr,lpart)
      common /iparts/ ipart(li,lpart)
      common /iparti/ n,nr,ni
      logical ifpts
      
      
      if(istep.eq.0) then
      nr=lr
      ni=li
      call rzero(rpart,lr*lpart)
      call izero(ipart,li*lpart)
      call particle_init       (rpart,nr,ipart,ni,n)
      write(6,*) n,'number of particles'
c      call particle_find       (rpart,nr,ipart,ni,n)
      else
      call particle_advect_up  (rpart,nr,ipart,ni,n)
      endif
      call concentration (rpart,nr,ipart,ni,n)
      write(6,*) 'made it here too'
c      call exitt
c      ipstep = iostep ! Particle output coordinated with iostep
c       ipstep = iostep/(iostep/10)
c       if(mod(istep,ipstep).eq.0.or.istep.eq.0)
      call particle_out        (rpart,nr,ipart,ni,n)
       
c      call exitt
c      if(istep.eq.3) call exitt
      return
      end
c-----------------------------------------------------------------------
      subroutine particle_init(rpart,nr,ipart,ni,n)

c     This version distributes particles through out the entire domain
c     at specified volume fraction 
c
      include 'SIZE'
      include 'TOTAL'
      
      real    rpart(nr,n)
      integer ipart(ni,n)
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      real wrk(lx1*ly1*lz1*lelt,ldim)
      common /outtmp/ wrk
      integer nw,xline,yline,zline
      real xp,zp,yp,pd,dx,dz,dy
      
      
      integer lcount,icount
      save    lcount,icount
      data    lcount,icount /0,0/
      
      call set_part_pointers
c     Particle volume fraction calculation
!//////////////////////////////////////////////////////////////////     
c      real vf,pvt,num_part,dx,dy,dz
c     real xp,yp,zp
      call domain_size(xmin,xmax,ymin,ymax,zmin,zmax) 
c      pi    = 4.*atan(1.0)
      pd = .01
c      vf=0.0005
c      pv=(4/3)*pi*(pd/2)**3
c      !particle volume
c      pvt=volvm1*vf
c      num_part=nint(pvt/pv)
c      write(6,2) pvt,volvm1,num_part
c 2    format(1p1e15.3,'particle volume', 1p1e15.3, 'total volume',
c     $ 1p1e15.3,'number of particles') 
!!////////////////////////////////////////////////////////////////////     
c      call exitt
      
c      if (mod(istep,ipstep).ne.0) return
c      write(6,*) nr,'nr'
c      write(6,*) ni, 'ni'
      llpart = lpart
     
      k  = icount       ! icount = total # particles emitted
      l  = lcount       ! Index into local particle array
      dx=0.2450
      dy=0.22
      dz=dx
      nw = 5
      xp=0.0
      yp=xp
      zp=yp
      do xline=1,nw
         xp =  xmin+pd+(xline-1)*dx
      	do yline = 1,10   ! 2 lines at different heights
        	yp = ymin+pd+(yline-1)*dy
        	do zline = 1,nw  ! nw points on a wire
         		zp = zmin+pd + (zline-1)*dz

         if (mod(k,np).eq.nid) then ! Inject particle _only_ on this processor
            l=l+1  ! local count
            rpart(jx,l) = xp
            rpart(jy,l) = yp
            rpart(jz,l) = zp
            rpart(jf0,l)= 0.0
            rpart(jar,l)= 0.0
            ipart(jai,l)=k+1
         endif
         k = k+1     ! Total count
        	enddo
     	enddo
      enddo
      write(6,*) k ,'total particle count'
      icount = k
      lcount = l
      npart=0
      npart  = max(n,lcount)
      n= npart
      write(6,*) n, 'number of particles per procs'
      call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),vx,vy,vz)
      nxyz  = nx1*ny1*nz1
      ntot  = nxyz*nelt
c      write(6,'(1pe15.7)') (wrk(n,3),k=1,ntot)
c      call exitt
      call particle_interp (rpart,nr,ipart,ni,n,wrk)
c      write(6,*) 'Didnt make it here'
      ! initialize particles velocities
      !do i=1,n
      !rpart(jv0+1,n)=-1.0
      !enddo
      return
      end
c---------------------------------------------------------------------------
      subroutine collisions_detection(rpart,nr,n)
      include 'SIZE'
      include 'TOTAL'
      
      real    rpart(nr,n)
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      common /maxmin /xyzmax(ldim),xyzmin(ldim)
      integer e,f,eg,m,icalld,post
      real yar(lx1*ly1*lz1*lelt),zar(lx1*ly1*lz1*lelt),
     &     xar(lx1*ly1*lz1*lelt),xyzmax,xyzmin,pts,vel
      real normc,wallpt,tolin,vold
      real incr(ldim,lpart)
      save icalld
      data icalld /0/
      logical log
      
      if (icalld.eq.0) then
      j=0
      nface = 2*ndim
      do e=1,nelt
      do f=1,nface
         if ((cbc(f,e,1).eq.'W  '))  then
          eg=lglel(e)
          call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
             
             do iz=kz1,kz2
             do iy=ky1,ky2
             do ix=kx1,kx2
             j=j+1
             xar(j)=xm1(ix,iy,iz,e)
             yar(j)=ym1(ix,iy,iz,e)
             zar(j)=zm1(ix,iy,iz,e)
             enddo
             enddo
             enddo

         endif 
      enddo
      enddo
      xyzmax(1)=glmax(xar,j)
      xyzmin(1)=glmin(xar,j)
      xyzmax(2)=glmax(yar,j)
      xyzmin(2)=glmin(yar,j)
      xyzmax(3)=glmax(zar,j)
      xyzmin(3)=glmin(zar,j)
c      write(6,*)xyzmax(1),xyzmin(1),xyzmax(2),xyzmin(2),
c     & xyzmax(3),xyzmin(3) 
c      call exitt
      icalld=1
      endif
      
c      do n=1,npart
c         write(6,*), uu(2,n), 'v_collision1'
c      enddo   
      post=0.0  
      write(6,*) 'Checking for collisions...'
      jx0 = jx
      do i=1,n
      	do j=0,ndim-1
      	log=.false.
      		if(rpart(jx0+j,i).gt.xyzmax(j+1)) then 
      			normc=rpart(jx0+j,i)
      			wallpt=xyzmax(j+1)
          		rpart(jx0+j,i)=wallpt-abs(normc-wallpt) 
      			vold=rpart(jv0+j,i)
      			rpart(jv0+j,i)=-rpart(jv0+j,i)
      			rpart(jv1+j,i)=-rpart(jv1+j,i)
      			rpart(jv2+j,i)=-rpart(jv2+j,i)
      			rpart(jv3+j,i)=-rpart(jv3+j,i)
      			if(j.eq.1) then
      			rpart(jx0+j,i)=xyzmin(j+1)+abs(normc-wallpt)
      			vold=rpart(jv0+j,i)
      			rpart(jv0+j,i)=-rpart(jv0+j,i)
      			rpart(jv1+j,i)=-rpart(jv1+j,i)
      			rpart(jv2+j,i)=-rpart(jv2+j,i)
      			rpart(jv3+j,i)=-rpart(jv3+j,i)
      			write(6,*) 'recycled pt'
      			endif 
      			log=.true.
      		elseif (rpart(jx0+j,i).lt.xyzmin(j+1)) then
      			normc=rpart(jx0+j,i)
      			wallpt=xyzmin(j+1)
      			rpart(jx0+j,i)=wallpt+ abs(normc-wallpt)
      			vold=rpart(jv0+j,i)
      			rpart(jv0+j,i)=-rpart(jv0+j,i)
      			rpart(jv1+j,i)=-rpart(jv1+j,i)
      			rpart(jv2+j,i)=-rpart(jv2+j,i)
      			rpart(jv3+j,i)=-rpart(jv3+j,i)
      			log=.true.
      		endif
     		
      	if(log.and.post.le.5) then
      write(6,'(a,e15.7,i5.1,e15.7,a)') 
     &     ' Collison detected the point normal to the wall is:',
     &     normc,j+1,wallpt,' Norm pt, dim, wall pt'
      write(6,'(a,e15.7,a,2e15.7)') 
     &     ' New norm pt is: ',
     &     rpart(jx0+j,i),' velocity change ',rpart(jv0+j,i),vold
        !else
         post=post+1
        !write(6,*) 'No collisions...'
        endif
c      vnorm=vel(i,2,n)
      !do some calls for the new positions and velocities
           enddo 
      enddo
c      do n=1,npart
c         write(6,*), uu(2,n), 'v_collision2'
c      enddo   
c       call exitt
      return
      end          
c-------------------------------------------------------------------------
      subroutine particle_advect_up(rpart,nr,ipart,ni,n)
      include 'SIZE'
      include 'TOTAL'
     
      real    rpart(nr,n)
      integer ipart(ni,n)
      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      real  wrk(lx1*ly1*lz1*lelt,ldim)  
      common /outtmp/ wrk
      real g(3),tau,c(3),rha   
      save g
      data g /9.80665,0.0,0.0/
      
      tau=0.001
      call particle_AB3_coef(c)
      jx0=jx
c      write(6,*) rpart(jv0+1,1), 'velocity'
c      write(6,*) rpart(jx0+1,1), 'y-position'
c           Solve for velocity at time t^n
      do j=0,ndim-1
      do i=1,n
         rpart(ju3+j,i)=rpart(ju2+j,i)
         rpart(ju2+j,i)=rpart(ju1+j,i)
         rpart(ju1+j,i)=rpart(ju0+j,i)
         rpart(jv3+j,i)=rpart(jv2+j,i)
         rpart(jv2+j,i)=rpart(jv1+j,i)
         rpart(jv1+j,i)=rpart(jv0+j,i)
         rpart(jx3+j,i)=rpart(jx2+j,i)
         rpart(jx2+j,i)=rpart(jx1+j,i)
         rpart(jx1+j,i)=rpart(jx0+j,i)
      enddo
      enddo
      
      ! Acceleration of particle
      do i=1,n
      do j=0,ldim-1     
      	rpart(ja0+j,i)=(rpart(ju1+j,i)-rpart(jv1+j,i))/tau-g(j+1)
        !rha=(-rpart(jv1+j,i))/tau-g(j+1)
c        if(istep.eq.0) rha=0.0
      	rpart(jv0+j,i) = rpart(jv1+j,i)+
     $               dt*(c(1)*rpart(ja0+j,i)
     $               +c(2)*rpart(ja1+j,i) 
     $               + c(3)*rpart(ja2+j,i))
      	rpart(jx0+j,i) =rpart(jx1+j,i)+
     $              dt*(c(1)*rpart(jv0+j,i)
     $ 			    + c(2)*rpart(jv1+j,i)
     $  			+ c(3)*rpart(jv2+j,i))
      	
      enddo	
      enddo
	  
	  call collisions_detection (rpart,nr,n)
	  call opcopy(wrk(1,1),wrk(1,2),wrk(1,3),vx,vy,vz)
      call particle_interp (rpart,nr,ipart,ni,n,wrk)
      
c      filename='vel_pos.txt'
c      acc='Append'
c      fr='formatted'
c      write(6,*) rpart(jv0+1,1), 'velocity'
c      write(6,*) rpart(jx0+1,1), 'y-position'
      
      
c      if(istep.eq.0) acc='Sequential'
c        write(6,*) filename, 'Is this the problem?'
c        open(unit=72,file=filename,form=fr,access=acc)
c             write(72,2) uu(2,1),x1(2,1),time
c 2             format ('',3e20.12)   
      
c      if(istep.eq.iostep) close(72)
      
      return 
      end       
c-----------------------------------------------------------------------
      subroutine particle_interp (rpart,nr,ipart,ni,n,wrk)
      include 'SIZE'
      include 'TOTAL'
      
      real    rpart(nr,n)
      integer ipart(ni,n)

      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      parameter (lrf=18*ldim+2,lif=6)
      real uvw(1),wrk(1)
      common /iv_intp/ ihandle
      
      
      call particle_find(rpart,nr,ipart,ni,n)
      
      ltot = lelt*lx1*ly1*lz1
      nflds  = ndim
c      write(6,*) n,'is this the right number of particles?'
      	do ifld = 1,nflds
         	iin    = (ifld-1)*ltot + 1
         	iout   = (ifld-1)*n + 1
         	is_out = 1
            iout   = ifld
           	is_out = nflds
c           	write(6,*) iin,'iin'
c           	write(6,*) wrk(iin+1) , 'vel' 
c           	write(6,*) iout,'iout' 
c          	write(6,*) is_out, 'is_out'
         	call findpts_eval(ihandle,
     &    	               rpart(ju0+iout-1,1),lrf,
     &                     ipart(jrc,1),lif,
     &                     ipart(jpt,1),lif,
     &                     ipart(je0,1),lif,
     &                     rpart(jr ,1),lrf,n,
     &                     wrk(iin))
      	enddo
c      write(6,*)n,'made I made it here?'
c      do i=1,n
c      do j=0,ldim-1
c      	rpart(ju0+j,i)=uvw(j+1,i)
c      	
c      enddo	
c      enddo
       if(istep.eq.0) then
       
       do i=1,n
       do j=0,ldim-1
       rpart(jv0+j,i)=rpart(ju0+j,i) 
       enddo
       enddo
       endif
c      do i=1,n
c      write(6,'(A,1p3e15.7,I2)')'V',(rpart(jv0+k,i),k=0,ndim-1)
c     & ,ipart(jrc,i)
c      write(6,'(A,1p3e15.7,I2)')'U',(rpart(ju0+k,i),k=0,ndim-1)
c     & ,ipart(jrc,i)
c      enddo
c      write(6,*) n, nr
c      write(6,*) rpart(ju0+1,1),uvw(2)
     
      write(6,*) 'particle interpolation complete'
c      call exitt
      return 
      end
c--------------------------------------------------------------------- 
      subroutine particle_AB3_coef(c)
	  !Adam Bashforth Coefficients 
      include 'SIZE'
      include 'TOTAL'
      real c(3)
      
      if (istep.le.1) then      ! AB1
         c(1) = 1.
         c(2) = 0.
         c(3) = 0.
      elseif (istep.eq.2) then  ! AB2
         c(1) = 3
         c(2) = -1.
         c(3) = 0
         c(1) = c(1)/2.
         c(2) = c(2)/2.
      else                    ! AB3
         c(1) = 23.
         c(2) = -16.
         c(3) = 5.
         c(1) = c(1)/12.
         c(2) = c(2)/12.
         c(3) = c(3)/12
      endif
      
      return
      end
c----------------------------------------------------------------------
      subroutine particle_out (rpart,nr,ipart,ni,n)
      include 'SIZE'
      include 'TOTAL'
      
      real    rpart(nr,n)
      integer ipart(ni,n)

      real x(ldim,lpart),partv(lpart)
      common /scrns/ x_tmp(ldim+2,lpart),work(ldim+2,lpart)
     $              ,v_tmp(ldim+1,lpart)
      character*128 fname

      common /ptpointers/ jrc,jpt,je0,jps,jai,nai,jr,jd,jx,jy,jz,jx1
     $               ,jx2,jx3,ja0,ja1,ja2,ja3,jv0,jv1,jv2,jv3
     $				 ,ju0,ju1,ju2,ju3,jf0,jar,nar
      integer icalld
      save    icalld
      data    icalld  /-1/

      icalld = icalld+1
      if (nid.eq.0) then
        write(fname,1) icalld
 1      format('part',i5.5,'.3D')
        open(unit=72,file=fname)

        write(fname,2) icalld
 2      format('vel',i5.5,'.3D')
        open(unit=73,file=fname)
      endif

      npt_total = iglsum(n,1)
      npass = npt_total / lpart
      if (npt_total.gt.npass*lpart) npass = npass+1
      ilast = 0
      do ipass=1,npass
        mpart = min(lpart,npt_total-ilast)
        i0    = ilast
        i1    = i0 + mpart
        ilast = i1

        call rzero(x_tmp,(ldim+2)*lpart)
        call rzero(v_tmp,(ldim+1)*lpart)
        do ii=1,n ! loop over all particles
          if (i0.lt.ipart(jai,ii).and.ipart(jai,ii).le.i1) then
            i = ipart(jai,ii)-i0
            call copy(x_tmp(1,i),rpart(jx,ii),ldim)  ! Coordinates
            x_tmp(ldim+1,i) = ipart(jpt,ii) ! MPI rank
            x_tmp(ldim+2,i) = ipart(jai,ii) ! Part id 
            call copy(v_tmp(1,i),rpart(jv0,ii),ldim)  ! Velocity 
            v_tmp(ldim+1,i) = ipart(jai,ii) ! Part id 
          endif
        enddo
        call gop(x_tmp,work,'+  ',(ldim+2)*lpart)
        call gop(v_tmp,work,'+  ',(ldim+1)*lpart)
        if (nio.eq.0) write(72,3)((x_tmp(k,i),k=1,ldim+2),i=1,mpart)
        if (nio.eq.0) write(73,4)((v_tmp(k,i),k=1,ldim+1),i=1,mpart)
 3      format(1p5e17.9)
 4      format(1p4e17.9)
      enddo

      if (nio.eq.0) close(72)
      if (nio.eq.0) close(73)
      return
      end
C-----------------------------------------------------------------------c
c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)
      return
      end
c
c automatically added by makenek
      subroutine cmt_switch ! to set IFCMT logical flag
      include 'SIZE'
      include 'INPUT'
      IFCMT=.false.
      return
      end
c
c automatically added by makenek
      subroutine usrflt(rmult) ! user defined filter
      include 'SIZE'
      real rmult(lx1)
      call rone(rmult,lx1)
      return
      end
c
c automatically added by makenek
      subroutine userflux ! user defined flux
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      real fluxout(lx1*lz1)
      return
      end
c
c automatically added by makenek
      subroutine userEOS ! user defined EOS 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
