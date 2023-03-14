      program lisalanc
      implicit none
      integer nss,nsp11,nblockk,ntott,nmpara,nmaxx,Iwmax
      integer Iwmaxreal,blocksize1,prozessoren
      integer ns,nsp1,nblock,ntot,nso2,i,j,k,nmax,nchoos
      integer imaxmu,iteramax,imu,itera,Nitermax,ksteps,ikx,iky
      logical normOrNot
!Parameters for magnetic field
      integer L   !size of magnetic unit cell
      real*8 Bm  !magnetic field B, should be dreal(p/L)
!For test one can set B=0 with finite L and should obtain the normal DMFT results
      integer omnumber,omrest,omoffset
      real*8 deltamu,testdmft,piob,xmu0,chi2,chi22,kx,ky, alpha
c      real*8 time,mclock
      real*8 zpart
c        in the whole file the only changes are: init.h changed place, number of processors and number of frequencies!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      integer*4 unlink,a
      include'init.h'                                       !include.init has to stand here so that the code does not give NaN
      !parameter (nss=5)
      parameter (prozessoren=36)                            !changed it from 400 to 36
      parameter (nsp11 =nss+1) 
      parameter (nblockk=nsp11*nsp11)
      parameter (ntott=4**nss)
      parameter (nmpara=30)                                                                                                                                    
      parameter (ksteps=100)  
      parameter (alpha=0.1) 
      parameter (normOrNot= .True.)                                                                                                                                
!Magnetic field parameter definition
      parameter (L=3)
      parameter (Bm=dfloat(1)/dfloat(L))
      !include'init.h'                                        !This had to be commented out, since used above
      parameter (Iwmax= 32768)                                  !Use 2**5 to be faster for testing
      parameter (Iwmaxreal =2048)
      real*8 tpar(nss),epsk(nss),uhub,hmag,xmu,beta,densimp,pi
      real*8 wlow,wup,deltino,range,sumdos,help1,help2,dens,sum
      real*8 om(0:Iwmax), dos(0:Iwmaxreal),dospart(0:Iwmaxreal)
      real*8 omm
      complex*16 ni(-10*Iwmax:10*Iwmax)
      real*8 np(0:Iwmax)
      integer offsetnum(0:prozessoren-1),recvnumber(0:prozessoren-1)
      integer offsetnuml(0:prozessoren-1),recvnumberl(0:prozessoren-1)
      integer offsetnum1(0:prozessoren-1),recvnumber1(0:prozessoren-1)
      integer offsetnume(0:prozessoren-1),recvnumbere(0:prozessoren-1)
      integer lehmnumber,eigennumber
      
      complex*16 omr(0:Iwmaxreal)
      integer ireal
      complex*16 Xi,c0,c1,ci
      complex*16 Gw(0:Iwmax),self(0:Iwmax),Gwpart(0:Iwmax)
      complex*16 Gww(0:2*Iwmax-1)
      complex*16 Gwreal(0:Iwmaxreal),Gwrealpart(0:Iwmaxreal)
      complex*16 sigre(0:Iwmaxreal)
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax),G0wwand(0:2*Iwmax-1)
      complex*16 G0wpart(0:Iwmax)
!Orbital matrices for magnetic field calc
      complex*16 G0wmat(0:Iwmax,1:L,1:L),G0wpartmat(0:Iwmax,1:L,1:L)
      complex*16 diff(1:L)
c     matrix & friends
      integer ioffset(nblockk+1)
      integer b(nss)
      integer ntable(ntott)
      integer*1 abas(nss)
      integer nleng(nblockk)
      integer*1 imat(nmaxx,nmaxx,nblockk)
      integer idmat(ntott)
      integer iddomat(ntott)
      integer*1 idouble(ntott)
c     for minimization
      real*8 hess(nmpara**2),g(nmpara),xtemp(nmpara),
     &     w(nmpara**2),xprmt(nmpara)
c     for dspev
      real*8 work(3*nmaxx)
      integer info
      real*8 eig(nmaxx),zeig(nmaxx,nmaxx)
      real*8 eigtotal(nmaxx,nblockk+1)
      real*8 eigold(nmaxx),zold(nmaxx,nmaxx)
      real*8 realmat(nmaxx*(nmaxx+1)/2)
      real*8 xmat(nmaxx,nmaxx)
      real*8 rlehm(nmaxx,nmaxx)
      real*8 rlehmtotal(nmaxx,nmaxx,nblockk+1)
      character*80 xyz 

      complex*16 omega, adummy, cdummy1

c      time=mclock()/1
      real*8 Epskold(nss),tparold(nss)

      complex*16 sig0
      real*8 zed,dens_wanted,diffdens,chi,densold,emin,dostest,emin1
      real*8 xmularge,xmusmall,denslarge,denssmall,densmix,emintot
      real*8 doubletot,zparttot
      integer ifix,inew,iauto
      integer iattempt,idritto,ilarge,ismall
      integer imix,iteraok,iexp
      real*8 threshold,th0,b1,b2,E
      real*8 doublep(ntott)
      real*8 double
      
      include 'mpif.h'
      integer ierror,myid,nprocs,mpi_type,myblock
      integer sendstatus(MPI_STATUS_SIZE)
      integer recvstatus(MPI_STATUS_SIZE)
      integer request,next,previous
      real*8 Energy
      external Energy
      real*8 Energya(0:ksteps,0:ksteps)

      


c      complex*16 gamma(-Iwmax:Iwmax,-Iwmax:Iwmax,-Iwmax:Iwmax)
c      complex*16 gammach(-Iwmax:Iwmax,-Iwmax:Iwmax,-Iwmax:Iwmax)
c      complex*16 Chi0inv(-Kmax:Kmax-1,-Kmax:Kmax-1,
c     $     -Iwmax:Iwmax-1,-Iwmax:Iwmax-1)


      logical bethe,twodim,symm,gwcalc
     
      common /symm/ symm  
       
      
      call MPI_INIT(ierror)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierror)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
  
      myblock=myid+1
      next=myid+1
      previous=myid-1
      
c      if (myid.eq.(nprocs-1)) then
c         next=0
c      endif

c      if (myid.eq.0) then
c         previous=nprocs-1
c      endif

      omnumber=(Iwmax+1)/nprocs
      omrest=mod(Iwmax+1,nprocs)
      do i=0,nprocs-1
         if (i.lt.omrest) then
            recvnumber(i)=omnumber+1
            offsetnum(i)=i*(omnumber+1)
         else
            recvnumber(i)=omnumber
            offsetnum(i)=omrest*(omnumber+1)+(i-omrest)*omnumber
         endif
      enddo
      do i=0,nprocs-1
         recvnumber1(i)=recvnumber(i)*L**2
      enddo
      offsetnum1(0)=0
      do i=1,nprocs-1
         offsetnum1(i)=offsetnum1(i-1)+recvnumber1(i-1)
      enddo
      if (myid.lt.omrest) then
         omnumber=omnumber+1
      endif
      if (myid.lt.omrest) then
         omoffset=myid*omnumber-1
      else
         omoffset=omrest*(omnumber+1)+(myid-omrest)*omnumber-1
      endif
      
      
      bethe=.false.
      twodim=.true.
      symm=.false.

      if (myid.eq.0) then
         if (bethe) then
            write(6,*) 'BETHE lattice case'
         else
            if (twodim) then
               write(6,*) '2D-simple cubic lattice'
            else
               write(6,*) '3D cubic lattice'
            endif
         endif
      endif 

      Pi=dacos(-1.d0)
      Xi=dcmplx(0.d0,1.d0)

      do ikx=0,ksteps
         kx=Pi*dfloat(ikx)/dfloat(ksteps)
         do iky=0,ksteps
            ky=Pi*dfloat(iky)/dfloat(ksteps)
            Energya(ikx,iky)=Energy(kx,ky)
         enddo
      enddo

      if (myid.eq.0) then
         open (30,file='hubb.dat',form='formatted',status='old')
         open (15,file='hubb.andpar',form='formatted',status='old')
         
         call datain(nss,uhub,xmu,hmag,ns,beta,wlow,wup,deltino,
     $               imaxmu,deltamu,iteramax,testdmft,dens_wanted,ifix,
     $               inew,iexp,th0,iauto)
      endif
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call flush(6)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
c      call MPI_BCAST(nss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(uhub,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(hmag,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(ns,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(beta,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(wlow,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(wup,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(deltino,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(imaxmu,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(deltamu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(iteramax,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(testdmft,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(dens_wanted,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(ifix,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(inew,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(iexp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(th0,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call MPI_BCAST(iauto,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierror)

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      if (myid.eq.((ns+1)**2-1)) then
         next=0
      endif

      if (myid.eq.0) then
         previous=(ns+1)**2-1
      endif

      if (myblock.gt.(ns+1)**2) then
         myblock=1
      endif

      piob=Pi/beta
      nsp1=ns+1
      nblock=nsp1*nsp1
      ntot=4**ns
      nso2=ns/2
      diffdens=10000.0d0
      double=0.0d0
      zpart=0.0d0
      dens=0.0d0
      emin1=10000.0d0

      if (myid.le.nblock-1) then
         lehmnumber=nmaxx**2
         eigennumber=nmaxx
      else 
         lehmnumber=1
         eigennumber=1
      endif

      do i=0,nprocs-1
         if (i.le.nblock-1) then
            recvnumberl(i)=nmaxx**2
            offsetnuml(i)=i*nmaxx**2
            recvnumbere(i)=nmaxx
            offsetnume(i)=i*nmaxx
         else
            recvnumberl(i)=1
            offsetnuml(i)=nblock*nmaxx**2+i-nblock
            recvnumbere(i)=nmaxx
            offsetnume(i)=nblock*nmaxx+i-nblock
         endif
      enddo

      do i=1,nmaxx
         do j=1,nmaxx
            rlehm(i,j)=0.0d0
         enddo
      enddo

      do i=1,nmaxx
         eig(i)=0.0d0
         do j=1,nmaxx
            zeig(i,j)=0.0d0
         enddo
      enddo

      b(1)=1
      do i=2,ns
         b(i)=b(i-1)*4
      enddo

      do i=1,nss
         tpar(i)=0.d0
         epsk(i)=0.d0
      enddo

      range=(wup-wlow)/dfloat(Iwmaxreal)
      do i=0,Iwmaxreal
         omr(i)=wlow+range*dfloat(i)+Xi*deltino
      enddo

      nmax=nchoos(ns,nso2)*nchoos(ns,nso2)

      if (myid.eq.0) then
         write(6,*) 'nmax', nmax, 'nmaxx', nmaxx
         if (nmax.ne.nmaxx) then
            write(6,*)'You are running serious troubles'
            write(6,*)'put nmaxx =',nmax
            stop
         endif
         write(6,*)'U  = ',uhub
         write(6,*)'mu = ',xmu
         write(6,*)'h = ',hmag
         write(6,*)'beta =',beta
         write(6,*) 'Iwmax =',Iwmax
         write(6,*)'for real frequencies: ',wlow,wup
         call flush(6) 
         if (ifix.eq.0) then
            write(6,*)'Working at fixed mu =',xmu
         elseif(ifix.eq.1) then
            write(6,*)'Working at fixed density =',dens_wanted  
            if(beta.ge.200.d0) then
               chi= 1.2d0/(1+0.6d0*dabs(uhub))
               xmu=(dens_wanted-1.d0)/chi
               xmu=0.5d0*uhub+xmu
            endif
         else
            write(6,*)'ifix should be 0, 1 !!'
            stop
         endif
         chi= 1.2d0/(1+0.6d0*dabs(uhub))
         call flush(6)
         write(6,*)'Initial mu =',xmu
         open(34,file='self-en_wim',form='formatted',status='unknown')
         open(90,file='gm_wim',form='formatted',status='unknown')
         open(91,file='gm_wre',form='formatted',status='unknown')
         open(92,file='g0m',form='formatted',status='unknown')
         open(93,file='g0mand',form='formatted',status='unknown')
         open(99,file='gw_cut',form='formatted',status='unknown')
         open(79,file='vert_chi',form='formatted',status='unknown')
         open(35,file='vsmu',form='formatted',status='unknown')
         open(89,file='GAMMA_WIM',form='formatted',status='unknown')
         open(85,file='SELF_LOC',form='formatted',status='unknown')
         open(50,file='SELF_Q',form='formatted',status='unknown')
         open(51,file='chiPI_wim',form='formatted',status='unknown')
         open(52,file='chiPiPi',form='formatted',status='unknown')
         open(74,file='chiSC_test',form='formatted',status='unknown')
         open(73,file='chiSC_wim',form='formatted',status='unknown')
         open(72,file='besetungsgsz',form='formatted',status='unknown')
         open(71,file='zustandsdicht',form='formatted',status='unknown')
         open(43,file='testdensity',form='formatted',status='unknown')
         open(47,file='testg0',form='formatted',status='unknown')
         open(78)
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
      call MPI_BCAST(chi,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
      call computematrix(ioffset,nblock,ntot,b,ns,abas,nsp1,nmax,
     $                   nleng,imat,idmat,iddomat,ntable,idouble)

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      
      if (myblock.eq.1) then
         blocksize1=1
      else
         blocksize1=nleng(myblock-1)
      endif

      if (myid.eq.0) then
         write(6,*)'# of independent blocks:',nblock
      endif

      do i=0,Iwmax
         om(i)=(2.d0*dfloat(i)+1.d0)*piob
      enddo
      
      do i=-10*Iwmax,10*Iwmax
         ni(i)=dfloat(i)*piob*Xi
      enddo

      do i=0,Iwmax
         np(i)=(2.d0*dfloat(i)+1.d0)*piob
      enddo

      if (myid.eq.0) then
         call initial(epsk,tpar,ns,xmu)
         write(6,*)'------------------------------------------------'
         write(6,*) 'starting Anderson parameters '
         write(6,*)'------------------------------------------------'
         write(6,*) 'Eps(k) '
         do i=2,ns
            write(6,'(2f27.16)')epsk(i)
         enddo
         write(6,*)'V(k) '
         do i=1,ns-1
            write(6,'(2f27.16)')tpar(i)
         enddo
         if (ifix.eq.1.and.inew.eq.0) read(15,*)xmu
         if (inew.eq.0.and.iauto.eq.0) read(15,*)xmu
         write(6,'(f27.16,"   #chemical potential")')xmu
         do i=2,ns
            write(86,*)epsk(i)
         enddo
      endif

      iteraok=1
      gwcalc=.false.
      
      do itera=1,iteramax
         if (myid.eq.0) then
            imix=0
            threshold=th0/dfloat(iteraok)**iexp
            if (threshold.lt.1.d-6) threshold=1.d-6
            epsk(1)=-xmu
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(tpar,nss,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(epsk,nss,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         CALL MPI_BARRIER(MPI_COMM_WORLD,ierror)

         do i=0,Iwmax
            call calcg0(om(i),cdummy1,tpar,epsk,ns,Xi,Iwmax)
            g0wand(i)=cdummy1
         end do
         if (myid.eq.0) then   
            write(6,*)'------------------------------------------------'
            write(6,*)'   Iteration : ',itera
            write(6,*)'------------------------------------------------'
            write(6,*)'Threshold for n=',threshold
            if (ifix.ne.0) write(6,*)'initial mu = ',xmu
            call flush(6)
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(threshold,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         iattempt=0 
         idritto=0
         ilarge=0
         ismall=0
         double=0.d0

         if (myid.lt.(ns+1)**2) then

         call diag(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,          !Integer-Input   
     $              nleng(myblock),
     $              b,ntable,ioffset,imat(1,1,myblock),              !Integer-Array-Input
     $              xmu,uhub,hmag,                                   !Real-Input
     $              tpar,epsk,                                       !Real-Array-Input
     $              emin1,                                            !Real-Output
     $              eig,zeig)                                        !Real-Array-Output
  
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN,
     $                   0,MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            emin=emintot
            write(6,*)emin
            call flush(6)
         endif

         call MPI_BCAST(emin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         do i=1,nleng(myblock)
            eig(i)=eig(i)-emin
         enddo

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         if (myid.lt.(ns+1)**2) then

         call MPI_ISSEND(eig,nmaxx,MPI_REAL8,next,
     $                 98,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(eigold,nmaxx,MPI_REAL8,previous,
     $                 98,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)
         call MPI_ISSEND(zeig,nmaxx**2,MPI_REAL8,next,
     $                 99,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(zold,nmaxx**2,MPI_REAL8,previous,
     $                 99,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)

         call lehmann(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,          !Integer-Input   
     $                nleng(myblock),blocksize1,
     $                ioffset,idouble,idmat,                           !Integer-Array-Input
     $                xmu,uhub,hmag,beta,                              !Real-Input
     $                eig,zeig,eigold,zold,                            !Real-Array-Input
     $                double,zpart,dens,                               !Real-Output
     $                rlehm)                                           !Real-Array-Output
         
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            doubletot=doubletot/zparttot
            densimp=densimp/zparttot
            if (ifix.eq.1) then
               diffdens=densimp-dens_wanted
            endif
         endif

         if (ifix.eq.1) then
            call MPI_BCAST(diffdens,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (ifix.eq.0) then

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Computation of the Green-Function on the imaginary axis           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

            if (myid.eq.0) then
               write(6,*)'Green-Function started'
               call flush(6)
            endif

            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8,
     $                          rlehmtotal,recvnumberl,offsetnuml,
     $                          MPI_REAL8,MPI_COMM_WORLD,ierror)   
            call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8,
     $                          eigtotal,recvnumbere,offsetnume,
     $                          MPI_REAL8,MPI_COMM_WORLD,ierror)   
            call computegimag(omoffset,omnumber,                 !Integer-Input
     $                        nsp1,nmaxx,Iwmax,prozessoren,        
     $                        nleng,                             !Integer-Array-Input  
     $                        beta,                              !Real-Input
     $                        eigtotal,rlehmtotal,om,            !Real-Array-Input
     $                        Xi,                                !Complex-Input
     $                        Gwpart)                            !Complex-Array-Output 
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)
            call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber),
     $                       omnumber,MPI_DOUBLE_COMPLEX,Gw,recvnumber,
     $                       offsetnum,
     $                       MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,ierror)
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)

            if (myid.eq.0) then
               do i=0,Iwmax
                  Gw(i)=Gw(i)/zparttot
               enddo
               write(6,*)'Green-Function finished'
               call flush(6)
            endif
            call MPI_BARRIER(MPI_COMM_WORLD,ierror)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   End Green-Function Calculation                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         endif
         
         if ((ifix.eq.0).and.(myid.eq.0)) then
            write(6,*)'densimp = ',densimp
         endif
 
         if (ifix.eq.1) then 
 67         continue
            if (myid.eq.0) then
               write(6,*)'<n>=',densimp,'  <n>-n0=',diffdens
            endif
            if (dabs(diffdens).lt.threshold) then               
               if (myid.eq.0) then
                  write(6,*)'------------------------'
                  write(6,*)'converged at step ',itera
                  write(6,*)'------------------------'
                  call flush(6)
               endif
c           converged now
               if (myid.eq.0) then
                  write(6,*) 'Computing G_imag'
                  call flush(6)
               endif

               

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Computation of the Green-Function on the imaginary axis           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               if (myid.eq.0) then
                  write(6,*)'Green-Function started'
                  call flush(6)
               endif

               call MPI_BARRIER(MPI_COMM_WORLD,ierror)
               call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8,
     $                             rlehmtotal,recvnumberl,offsetnuml,
     $                             MPI_REAL8,MPI_COMM_WORLD,ierror)   
               call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8,
     $                             eigtotal,recvnumbere,offsetnume,
     $                             MPI_REAL8,MPI_COMM_WORLD,ierror)

               call computegimag(omoffset,omnumber,              !Integer-Input
     $                           nsp1,nmaxx,Iwmax,prozessoren,        
     $                           nleng,                          !Integer-Array-Input  
     $                           beta,                           !Real-Input
     $                           eigtotal,rlehmtotal,om,         !Real-Array-Input
     $                           Xi,                             !Complex-Input
     $                           Gwpart)                         !Complex-Array-Output 
               call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
               call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber),
     $                          omnumber,MPI_DOUBLE_COMPLEX,Gw,
     $                          recvnumber,offsetnum,
     $                          MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,
     $                          ierror)
               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

               if (myid.eq.0) then
                  do i=0,Iwmax
                     Gw(i)=Gw(i)/zparttot
                  enddo
                  write(6,*)'Green-Function finished'
                  call flush(6)
               endif
               call MPI_BARRIER(MPI_COMM_WORLD,ierror)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   End Green-Function Calculation                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

               iteraok=iteraok+1
               if (myid.eq.0) then
                  open(101)
                  do i=0,Iwmax
                     write(101,*)i,om(i),Gw(i)
                  enddo
               endif

            else
c           not converged now
               if (myid.eq.0) then
                  write(6,*)'-Not yet converged-'
                  call flush(6)
               endif
               if (diffdens.gt.0.d0) then
                  ilarge=1
                  xmularge=xmu
                  denslarge=densimp
               endif
               if (diffdens.lt.0.d0) then
                  ismall=1
                  xmusmall=xmu
                  denssmall=densimp
               endif
               if (ilarge*ismall.eq.0) then
                  if (myid.eq.0) then
                     write(6,*)'Still lacking a bracket'
                  endif
                  xmu=xmu-chi*diffdens
                  if (myid.eq.0) then
                     write(6,*)'Try with xmu=',xmu
                  endif
                  idritto=idritto+1
                  if (idritto.eq.6) then
                     chi=chi*2.d0
                  endif

                  if (myid.lt.(ns+1)**2) then
                  
                  call diag(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock, !Integer-Input   
     $                       nleng(myblock),
     $                       b,ntable,ioffset,imat(1,1,myblock),     !Integer-Array-Input
     $                       xmu,uhub,hmag,                          !Real-Input
     $                       tpar,epsk,                              !Real-Array-Input
     $                       emin1,                                   !Real-Output
     $                       eig,zeig)                               !Real-Array-Output
  
                  endif

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)


                  call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN,
     $                            0,MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     emin=emintot
                     write(6,*)emin
                     call flush(6)
                  endif

                  call MPI_BCAST(emin,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  do i=1,nleng(myblock)
                     eig(i)=eig(i)-emin
                  enddo

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  if (myid.lt.(ns+1)**2) then

                  call MPI_ISSEND(eig,nmaxx,MPI_REAL8,next,
     $                          98,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(eigold,nmaxx,MPI_REAL8,previous,
     $                          98,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
                  call MPI_ISSEND(zeig,nmaxx**2,MPI_REAL8,next,
     $                          99,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(zold,nmaxx**2,MPI_REAL8,previous,
     $                          99,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
         
                  call lehmann(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock, !Integer-Input   
     $                         nleng(myblock),blocksize1,
     $                         ioffset,idouble,idmat,                  !Integer-Array-Input
     $                         xmu,uhub,hmag,beta,                     !Real-Input
     $                         eig,zeig,eigold,zold,                   !Real-Array-Input
     $                         double,zpart,dens,                      !Real-Output
     $                         rlehm)                                  !Real-Array-Output
         
                  endif

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)

                  if (myid.eq.0) then
                     doubletot=doubletot/zparttot
                     densimp=densimp/zparttot
                     diffdens=densimp-dens_wanted
                  endif
                  
                  call MPI_BCAST(diffdens,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  go to 67
               else
                  call MPI_BCAST(denssmall,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  call MPI_BCAST(denslarge,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     write(6,*)'Interpolation between:'
                     write(6,*)'mu, n (1):',xmusmall,denssmall
                     write(6,*)'mu, n (2):',xmularge,denslarge
                  endif
                  xmu=xmusmall+(xmularge-xmusmall)*
     $                 (denssmall-dens_wanted)/
     $                 (denssmall-denslarge)
                  if (myid.eq.0) then
                     write(6,*)'interpolated xmu',xmu
                  endif

                  if (myid.lt.(ns+1)**2) then

                  call diag(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock, !Integer-Input   
     $                       nleng(myblock),
     $                       b,ntable,ioffset,imat(1,1,myblock),     !Integer-Array-Input
     $                       xmu,uhub,hmag,                          !Real-Input
     $                       tpar,epsk,                              !Real-Array-Input
     $                       emin1,                                   !Real-Output
     $                       eig,zeig)                               !Real-Array-Output
  
                  endif

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN,
     $                            0,MPI_COMM_WORLD,ierror)
                  if (myid.eq.0) then
                     emin=emintot
                     write(6,*)emin
                     call flush(6)
                  endif

                  call MPI_BCAST(emin,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  do i=1,nleng(myblock)
                     eig(i)=eig(i)-emin
                  enddo

                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)

                  if (myid.lt.(ns+1)**2) then

                  call MPI_ISSEND(eig,nmaxx,MPI_REAL8,next,
     $                          98,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(eigold,nmaxx,MPI_REAL8,previous,
     $                          98,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)
                  call MPI_ISSEND(zeig,nmaxx**2,MPI_REAL8,next,
     $                          99,MPI_COMM_WORLD,request,ierror)
                  call MPI_RECV(zold,nmaxx**2,MPI_REAL8,previous,
     $                          99,MPI_COMM_WORLD,recvstatus,ierror)
                  call MPI_WAIT(request,sendstatus,ierror)

                  call lehmann(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock, !Integer-Input   
     $                         nleng(myblock),blocksize1,
     $                         ioffset,idouble,idmat,                  !Integer-Array-Input
     $                         xmu,uhub,hmag,beta,                     !Real-Input
     $                         eig,zeig,eigold,zold,                   !Real-Array-Input
     $                         double,zpart,dens,                      !Real-Output
     $                         rlehm)                                  !Real-Array-Output

                  endif
         
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                  
                  call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)
                  call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM,
     $                 0,MPI_COMM_WORLD,ierror)

                  if (myid.eq.0) then
                     doubletot=doubletot/zparttot
                     densimp=densimp/zparttot
                     diffdens=densimp-dens_wanted
                  endif
                  
                  call MPI_BCAST(diffdens,1,MPI_REAL8,0,
     $                           MPI_COMM_WORLD,ierror)
                  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
                          
                  goto 67
               endif
            endif
         endif


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            rewind(78)
            do i=0,Iwmax
               write(78,*)om(i),Gw(i)
            enddo
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(Gw,Iwmax+1,MPI_DOUBLE_COMPLEX,0,
     $                  MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
 68      if (myid.eq.0) then
            write(6,*) 'SELF_CONSISTENCY LOOP'
            write(6,*) 'for the 2D-cubic' 
         endif
c            if(bethe) then
c               write(6,*) 'Bethe lattice' 
c               do i=0,Iwmax
c                  adummy=om(i)*Xi-Epsk(1)-Gw(i)/4.d0
c                  if (iauto.eq.0)then
c                     G0w(i)=1.d0/adummy
c                  else
c                     G0w(i) = (0.5d0/adummy+0.5d0*G0wand(i))
c                  endif
c               enddo 
c            else
c               if(twodim) then
c                  write(6,*) 'for the 2D-cubic' 

c         call selfconst2(omnumber,omoffset,Iwmax,ksteps,      !Integer-Input
c     $                   beta,xmu,                            !Real-Input  
c     $                   om,Energya,                          !Real-Array-Input
c     $                   Xi,                                  !Complex-Input      
c     $                   g0wand,Gw,                           !Complex-Array-Input
c     $                   G0wpart)                             !Complex-Array-Output

         call  selfconst2_magnetic(omnumber,omoffset,Iwmax,ksteps,L, !Integer-Input
     $                         beta,xmu,Bm,                          !Real-Input  
     $                         om,Energya,                           !Real-Array-Input
     $                         Xi,                                   !Complex-Input      
     $                         g0wand,Gw,                            !Complex-Array-Input
     $                         G0wpartmat)                           !Complex-Array-Output

!Here the local interacting (!) Green's function is written to G0w and G0wpart

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         do j=1,L
            do k=1,L
               call MPI_GATHERV(
     $              G0wpartmat(omoffset+1:omoffset+omnumber,j,k),
     $              omnumber,
     $              MPI_DOUBLE_COMPLEX,G0wmat(:,j,k),
     $              recvnumber,offsetnum,
     $              MPI_DOUBLE_COMPLEX,0,
     $              MPI_COMM_WORLD,ierror)
            enddo
         enddo
c               else  
c                  write(6,*) 'for the 3D-cubic'
c                  call selfconst3(g0wand,Gw,G0w,Xi,Pi,Iwmax,om,xmu)
c               endif
c            endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         if (myid.eq.0) then

            do j=1,L
               diff(j)=dcmplx(0.0d0,0.0d0)
            enddo
            
            do i=0,Iwmax
               do j=1,L
                  write(2000+j,*)om(i),
     $              dreal(G0wmat(i,j,j)),dimag(G0wmat(i,j,j))
               enddo
            enddo
            do i=0,Iwmax
               G0w(i)=G0wmat(i,1,1)
               do j=2,L
                  diff(j)=diff(j)+abs(G0w(i)-G0wmat(i,j,j))
               enddo
            enddo
            !write(*,*) sum(G0w)                                                                                            !TRy to find similarity between local and HLRN on 14.2.23
            !write(*,*) sum(G0wmat)
            do j=2,L
               write(6,*)"Difference to orbital 1 for orbital",j,diff(j)
            enddo

            do i=0,Iwmax
               G0w(i)=dcmplx(0.0d0,0.0d0)                                              !initialize G0w
            enddo

            do i=0,Iwmax
               do j=1,L
                  G0w(i)=G0w(i)+G0wmat(i,j,j)/dfloat(L)                                !average over all diagonale elements, which should be the same anyway
               enddo
            enddo

!Calculation of the new noninteracting GF of the AIM (also written to G0w)
            do i=0,Iwmax
               G0w(i)=1.0d0/G0w(i)+1.0d0/g0wand(i)-1.0d0/Gw(i)
!damping/mixing of self-consistency
               G0w(i)=(1-alpha)/G0w(i) + alpha*g0wand(i)                                !g0wand Greens of old one Anderson, G0w is new Greens (or other way around)     f_n = alpha f_(n-1) + (1-alpha) f_n
            enddo
         endif

         open(102)
         if (myid.eq.0) then
            do i=0,Iwmax
               write(102,*)om(i),dreal(G0w(i)),dimag(G0w(i))
            enddo
         endif

!         call MPI_FINALIZE(ierror)
!         stop

         if (myid.eq.0) then
            if (iauto.eq.0) then 
               write(6,*)'Average Double Occupancy :',doubletot
               goto 777
            endif
            
            Nitermax=4000
            write(6,*)'minimization'
            call flush(6)

            rewind(62)
            do i=0,Iwmax
               write(62,*) om(i),dreal(G0w(i)),dimag(G0w(i)),
     $                     dreal(Gw(i)),dimag(Gw(i))
               call flush(62)
            enddo
           
            call search(chi2,Nitermax,hess,g,xtemp,w,xprmt,nmpara,
     $           tpar,epsk,ns,piob,Xi,Iwmax,G0w,G0wand,om, normOrNot)
            call flush(6)
            rewind(15)

            rewind(63)
            do i=0,Iwmax
               write(63,'(3f17.10)') om(i),dreal(G0wand(i)),
     $              dimag(G0wand(i))
               call flush(63)
            enddo

            call rheader(15)
            
            read(15,*)

            do i=2,ns
               read(15,*)Epskold(i)
            end do

            read(15,*)

            do i=1,ns-1
               read(15,*)tparold(i)
            end do
            
            chi22=0.d0

            do i=1,ns-1
               chi22=chi22+(tpar(i)-tparold(i))**2
               chi22=chi22+(Epsk(i+1)-Epskold(i+1))**2
            enddo

            chi22=chi22/(2*ns-2)
            
            write(6,*)'------------------------------------------------'
            write(6,*) 'new Anderson parameters '
            write(6,*)'------------------------------------------------'
            write(6,*) 'Eps(k) '
            do i=2,ns
               write(6,'(2f27.16)')epsk(i)
            enddo
             write(6,*)'V(k) '
            do i=1,ns-1
               write(6,'(2f27.16)')tpar(i)
            enddo
            write(6,'(f27.16,"   #chemical potential")')xmu
            write(6,*)'------------------------------------------------'

            write(6,'(" convergence parameter :",e18.10)')chi22
            call flush(6)            

            write(22,*)itera,chi22
            call flush(22) 
            rewind(15)
            call wheader(15,ns,Iwmax)
            write(15,'(9a)')' Eps(k)'
            do i=2,ns
               write(15,*)Epsk(i)
            end do
            write(15,'(9a)')' tpar(k)'
            do i=1,ns-1
               write(15,*)tpar(i)
            end do
            write(15,*)xmu,'    #chemical potential'
            call flush(15)
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(chi22,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)

         if (chi22.lt.testdmft) goto 777
      enddo
C========+=========+=========+=========+=========+=========+=========+=$
c                               OUTPUT
C========+=========+=========+=========+=========+=========+=========+=$
         
777      continue

      if (iauto.ne.0) then
            
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_BCAST(tpar,nss,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(epsk,nss,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         call MPI_BCAST(xmu,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         
         if (myid.lt.(ns+1)**2) then

         call diag(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,          !Integer-Input   
     $              nleng(myblock),
     $              b,ntable,ioffset,imat(1,1,myblock),              !Integer-Array-Input
     $              xmu,uhub,hmag,                                   !Real-Input
     $              tpar,epsk,                                       !Real-Array-Input
     $              emin1,                                            !Real-Output
     $              eig,zeig)                                        !Real-Array-Output
  
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(emin1,emintot,1,MPI_REAL8,MPI_MIN,
     $                   0,MPI_COMM_WORLD,ierror)
         if (myid.eq.0) then
            emin=emintot
            write(6,*)emin
            call flush(6)
         endif

         call MPI_BCAST(emin,1,MPI_REAL8,0,MPI_COMM_WORLD,ierror)
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         do i=1,nleng(myblock)
            eig(i)=eig(i)-emin
         enddo

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

         if (myid.lt.(ns+1)**2) then

         call MPI_ISSEND(eig,nmaxx,MPI_REAL8,next,
     $                 98,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(eigold,nmaxx,MPI_REAL8,previous,
     $                 98,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)
         call MPI_ISSEND(zeig,nmaxx**2,MPI_REAL8,next,
     $                 99,MPI_COMM_WORLD,request,ierror)
         call MPI_RECV(zold,nmaxx**2,MPI_REAL8,previous,
     $                 99,MPI_COMM_WORLD,recvstatus,ierror)
         call MPI_WAIT(request,sendstatus,ierror)

         
         call lehmann(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,          !Integer-Input   
     $                nleng(myblock),blocksize1,
     $                ioffset,idouble,idmat,                           !Integer-Array-Input
     $                xmu,uhub,hmag,beta,                              !Real-Input
     $                eig,zeig,eigold,zold,                            !Real-Array-Input
     $                double,zpart,dens,                               !Real-Output
     $                rlehm)                                           !Real-Array-Output

         endif
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         call MPI_REDUCE(double,doubletot,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(zpart,zparttot,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)
         call MPI_REDUCE(dens,densimp,1,MPI_REAL8,MPI_SUM,
     $                   0,MPI_COMM_WORLD,ierror)

         if (myid.eq.0) then
            doubletot=doubletot/zparttot
            densimp=densimp/zparttot
         endif
         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Computation of the Green-Function on the imaginary axis           c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write(6,*)'Green-Function started'
            call flush(6)
         endif
         
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_ALLGATHERV(rlehm,lehmnumber,MPI_REAL8,
     $                       rlehmtotal,recvnumberl,offsetnuml,
     $                       MPI_REAL8,MPI_COMM_WORLD,ierror)   
         call MPI_ALLGATHERV(eig,eigennumber,MPI_REAL8,
     $                       eigtotal,recvnumbere,offsetnume,
     $                       MPI_REAL8,MPI_COMM_WORLD,ierror)   
         call computegimag(omoffset,omnumber,              !Integer-Input
     $                     nsp1,nmaxx,Iwmax,prozessoren,        
     $                     nleng,                          !Integer-Array-Input  
     $                     beta,                           !Real-Input
     $                     eigtotal,rlehmtotal,om,         !Real-Array-Input
     $                     Xi,                             !Complex-Input
     $                     Gwpart)                         !Complex-Array-Output 
         call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
         call MPI_GATHERV(Gwpart(omoffset+1:omoffset+omnumber),
     $                    omnumber,MPI_DOUBLE_COMPLEX,Gw,
     $                    recvnumber,offsetnum,
     $                    MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,
     $                    ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            do i=0,Iwmax
               Gw(i)=Gw(i)/zparttot
            enddo
            write(6,*)'Green-Function finished'
            call flush(6)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   End Green-Function Calculation                                    c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            do i=0,Iwmax-1
               Gww(i)=Gw(i)
               Gww(i+Iwmax)=1.d0/ni(2*(Iwmax+i)+1)
            enddo
           
            do i=0,2*Iwmax-1
               omm=dimag(ni(2*i+1))
               call calcg0(omm,cdummy1,tpar,epsk,ns,Xi,Iwmax)
               cdummy1=1.d0/cdummy1
               g0wwand(i)=cdummy1
            end do
         endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Computation of the Green-Function on the real axis                c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write(6,*)'Real Green-Function started'
            call flush(6)
         endif

         omnumber=(Iwmaxreal+1)/nprocs
         omrest=mod(Iwmaxreal+1,nprocs)
         do i=0,nprocs-1
            if (i.lt.omrest) then
               recvnumber(i)=omnumber+1
               offsetnum(i)=i*(omnumber+1)
            else
               recvnumber(i)=omnumber
               offsetnum(i)=omrest*(omnumber+1)+(i-omrest)*omnumber
            endif
         enddo
         if (myid.lt.omrest) then
            omnumber=omnumber+1
         endif
         if (myid.lt.omrest) then
            omoffset=myid*omnumber-1
         else
            omoffset=omrest*(omnumber+1)+(myid-omrest)*omnumber-1
         endif

         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call computegreal(omoffset,omnumber,                !Integer-Input
     $                     nsp1,nmaxx,Iwmaxreal,prozessoren,         
     $                     nleng,                             !Integer-Array-Input  
     $                     beta,                              !Real-Input
     $                     eigtotal,rlehmtotal,               !Real-Array-Input
     $                     omr,                               !Complex-Array-Input
     $                     Gwrealpart)                        !Complex-Array-Output

         call MPI_BARRIER(MPI_COMM_WORLD,ierror) 
         call MPI_GATHERV(Gwrealpart(omoffset+1:omoffset+omnumber),
     $                    omnumber,MPI_DOUBLE_COMPLEX,Gwreal,
     $                    recvnumber,offsetnum,
     $                    MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,
     $                    ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            do i=0,Iwmax
               Gwreal(i)=Gwreal(i)/zparttot
            enddo
            write(6,*)'Real Green-Function finished'
            call flush(6)
         endif
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   End Real-Green-Function Calculation                               c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         if (myid.eq.0) then
            write (6,*) 'Now computing Matsubara self-energy'       
            write(6,*) 'Pi is equal to ', Pi
         
            do i=0, Iwmax
               self(i) =g0wwand(i)-1.d0/Gw(i)
               write(34,'(3f17.10)') om(i),dreal(self(i)),dimag(self(i))
            enddo
         
            write (6,*) 'Now computing self-energy on the real axis'       
            
            do i=0,Iwmaxreal
               sigre(i)=(0.d0,0.d0)
               do j=1,ns-1
                  sigre(i)=sigre(i)-tpar(j)**2/(omr(i)-Epsk(1+j))
               end do
               write(22,'(4f17.10)') dreal(Gwreal(i)), dimag(Gwreal(i)),
     $              dreal(1.d0/Gwreal(i)),dimag(1.d0/Gwreal(i))
               sigre(i)=sigre(i)+omr(i)-Epsk(1)-
     $                  dconjg(Gwreal(i))/(Gwreal(i)*dconjg(Gwreal(i)))
            enddo
         endif
      endif

      if (myid.eq.0) then
         write(6,*)'Average Double Occupancy :',doubletot
      endif

      call MPI_BARRIER(MPI_COMM_WORLD,ierror)
      call MPI_BCAST(sigre,Iwmaxreal+1,MPI_DOUBLE_COMPLEX,0,
     $               MPI_COMM_WORLD,ierror)
      call MPI_BARRIER(MPI_COMM_WORLD,ierror)

      if(twodim) then
         open(57,file='dos.dat',form='formatted',status='unknown')
         call dos2(omnumber,omoffset,Iwmaxreal,ksteps,        !Integer-Input
     $             beta,xmu,Pi,                               !Real-Input  
     $             Energya,                                   !Real-Array-Input
     $             omr,sigre,                                 !Complex-Array-Input
     $             dospart)                                   !Real-Array-Output
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         call MPI_GATHERV(dospart(omoffset+1:omoffset+omnumber),
     $                    omnumber,MPI_REAL8,dos,
     $                    recvnumber,offsetnum,
     $                    MPI_REAL8,0,MPI_COMM_WORLD,
     $                    ierror)
         call MPI_BARRIER(MPI_COMM_WORLD,ierror)
         
         if (myid.eq.0) then
            sum=0.0d0
            do i=0,Iwmaxreal
               sum=sum+2.0d0*dos(i)*dreal(omr(1)-omr(0))/
     $            (dexp(beta*dreal(omr(i)))+1.0d0)
               write(57,'(3f17.10)') dreal(omr(i)), dos(i),sum
            enddo
            write(6,*) 'Final Spectral weight check=',sum
         endif
      endif
        
      if(bethe) then
         open(57,file='dos.dat',form='formatted',status='unknown')
         do i=0,Iwmaxreal
            dos(i)=0.d0
            do j=0,4999
               E=-1.d0+2.d0*dfloat(j)/dfloat(5000) +1.d0/dfloat(5000)
                  dos(i)=dos(i)-2.d0/(pi)**2*dsqrt(1.d0-E**2)*
     $                 dimag(1.d0/(omr(i)-Epsk(1)-E-sigre(i)))
     $                 /dfloat(5000)
            enddo
               if (myid.eq.0) then
                  write(57,'(2f17.10)') dreal(omr(i)),dos(i)
               endif
         enddo
      endif

      if (myid.eq.0) then
         write(6,*)'------------------------------------------------'
         write(6,*) 'final Anderson parameters '
         write(6,*)'------------------------------------------------'
         write(6,*)'cut from the following line...'
         write(6,*) 'Eps(k) '
         do i=2,ns
           write(6,'(2f27.16)')epsk(i)
         enddo
         write(6,*)'V(k) '
         do i=1,ns-1
            write(6,'(2f27.16)')tpar(i)
         enddo
         write(6,'(f27.16,"   #chemical potential")')xmu      
         write(6,*)'...to the line above and copy into hubb.andpar'
         write(6,*)'------------------------------------------------'
 
         do i=0,2*Iwmax-1
            write(90,'(3f27.16)') dimag(ni(2*i+1))
     $           ,dreal(Gww(i)),dimag(Gww(i))
         enddo
         
         do i=0,Iwmax
            G0w(i)=1.d0/G0w(i)
            write(92,'(3f17.10)') dimag(ni(2*i+1))
     $           ,dreal(G0w(i)),dimag(G0w(i))
         end do
               
         do i=0,2*Iwmax-1
            write(93,'(3f17.10)')dimag(ni(2*i+1)),dreal(G0wwand(i)),
     $           dimag(G0wwand(i))
         end do
         
         do i=0,Iwmaxreal
            write(91,'(5f17.10)')dreal(omr(i)),
     $           dreal(Gwreal(i)),dimag(Gwreal(i)),dreal(sigre(i)),
     $           dimag(sigre(i))
         end do
         
         sig0=g0w(0)-1.d0/gw(0)
         write(6,*) 'sig0 = ',sig0
         zed=1.d0-dimag(sig0)/piob
         zed=1.d0/zed
         write(6,*) 'zed = ', zed
         call flush(6)
          
         dostest=0.d0
         do i=1,ns-1
            dostest=dostest+tpar(i)**2
         enddo
         write(6,*) 'DOSTEST=', dostest

         write(35,888) uhub,xmu,densimp,emin,zpart
 888     format(f6.3,1x,f6.3,1x,f12.8,1x,f12.8,1x,f12.8,1x,f12.8)
         write(45,*) uhub,densimp, double
C========+=========+=========+=========+=========+=========+=========+=$
c     END OF Iteration 
C========+=========+=========+=========+=========+=========+=========+=$
         write(6,'(a20)')'     End lisalanc '
         write(6,'(a60)')'========================================'
      endif
      
      call MPI_FINALIZE(ierror)

      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine initial(epsk,tpar,ns,xmu)
      implicit none
      integer ns,i
      real*8 epsk(ns),tpar(ns),xmu

      rewind(15)
      call rheader(15)
c      Epsk(1)=-xmu
      read(15,*)
      do i=2,ns
         read(15,*)Epsk(i)
      end do
      read(15,*)
      do i=1,ns-1
         read(15,*)tpar(i)
      end do
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine rheader(k)
      implicit none
      integer i,k
 
      do 1 i=1,8
         read(k,*)
1     continue
      end 
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine wheader(k,ns,Iwmax)
      implicit none
      integer i,k
      integer ns,Iwmax

      character *80 xyz
      character *8 xxx
      xxx=' 1-band '
      write(k,'(a55)')'========================================'
      write(k,'(a25,a30)')xxx,'30-Sep-95 LANCZOS  '
      write(k,'(a55)')'========================================'
      rewind(30)
      write(k,'(4(a6,I4))') 'NSITE ',ns,'IWMAX',iwmax
      read(30,'(a60)')xyz
      read(30,'(a60)')xyz
      read(30,'(a60)')xyz
      do 3 i=1,4
      read(30,'(a60)')xyz
      write(k,'(a60)')xyz
3     continue 
      rewind(30)
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine search(fmin,Nitermax,hess,g,xtemp,w,xprmt,nmpara,
     $     tpar,epsk,ns,piob,Xi,Iwmax,G0w,G0wand,om, normOrNot)
      implicit none
      integer nmpara,ns,Iwmax,nbparm,icount,i,mode,Nitermax,iexit,nsym
      real*8 dfn,deps,fmin,hh
      real*8 hess(nmpara**2),g(nmpara),xtemp(nmpara),
     &     w(nmpara**2),xprmt(nmpara)
      real*8 tpar(ns),epsk(ns),piob  
      real*8 om(0:Iwmax)
      complex*16 Xi
      logical normOrNot
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      logical symm
      external difference
      data iexit/0/
      data hh/1.e-5/
      
      common /symm/ symm
 
c     number of parameters to be optimized (per spin):
c     eps(2)...... eps(ns) --->  ns-1
c     tpar(1)........ tpar(ns-1) --->  ns-1

      if(symm) then
         nsym=(ns-1)/2
         nbparm=2*(nsym-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
C
         icount=0
     
c         nsym=(ns-1)/2
         nbparm=2*(nsym)
         do i=2,nsym+1
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 11         continue
         end do
         
         do i=1,nsym
            icount=icount+1
            xtemp(icount)=tpar(i)
 12         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
      
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      else  
         nbparm=2*(ns-1)
         if (nbparm.gt.nmpara) stop ' nmpara too small'
C     
         icount=0
         do i=2,ns
            icount=icount+1	
            xtemp(icount)=Epsk(i)
 111        continue
         end do
         
         do i=1,ns-1
            icount=icount+1
            xtemp(icount)=tpar(i)
 112         continue
         end do
         
         do i=nbparm+1,nmpara
            xtemp(i)=0.d0
         end do
         
         do i=1,nbparm
            xprmt(i)=dabs(xtemp(i))+1.d-15
         end do
      endif

      

      mode=1
c     use va10 for search
c     make dfn and deps smaller to look if it helps, just an idea !!!!!!!!!!!!!!!!!!!!!
      dfn=-.5d0
      deps=.00001d0

      call minimize(difference,nbparm,xtemp,fmin,g,hess,w
     +     ,dfn,xprmt,hh,deps,mode,Nitermax,iexit,
     $     tpar,epsk,piob,Iwmax,Xi,ns,G0wand,G0w,om, normOrNot)
      write (6,30) iexit,fmin
c     do i=1,10 
c     write(73,*)xtemp(i),xprmt(i)
c     enddo
 30   format(' iexit fmin ',i5,e14.6)
      end
C========+=========+=========+=========+=========+=========+=========+=$
      subroutine minimize (funct, n, x, f, g, h, w, dfn, xm,
     $  hh, eps, mode, maxfn, iexit,tpar,epsk,piob,
     $     Iwmax,Xi,ns,G0wand,G0w,om, normOrNot)
      implicit none
      integer ns,Iwmax
      logical normOrNot
      integer np,n,n1,nn,is,iu,iv,ib,idiff,iexit,mode,ij,maxfn,i,j,
     $     i1,jk,ik,k,itn,ifn,link,int
      real*8 z,zz,dmin,f,df,dfn,aeps,eps,alpha,ff,tot,f1,f2,half,
     $     gys,dgs,sig,hh,gs0
      real*8 tpar(ns),epsk(ns),piob
      complex*16 Xi
      complex*16 G0w(0:Iwmax),G0wand(0:Iwmax)
      real*8  x(*), g(*), h(*), w(*), xm(*)
      real*8 om(0:Iwmax)
      external funct
      data half /0.5d0/


      np = n + 1
      n1 = n - 1
      nn=(n*np)/2
      is = n
      iu = n
      iv = n + n
      ib = iv + n
      idiff = 1
      iexit = 0
      if (mode .eq. 3) go to 15
      if (mode .eq. 2) go to 10
      ij = nn + 1
      do 5 i = 1, n
      do 6 j = 1, i
      ij = ij - 1
   6  h(ij) = 0.d0
   5  h(ij) = 1.d0
      go to 15
  10  continue
      ij = 1
      do 11 i = 2, n
      z = h(ij)
      if (z .le. 0.d0) return
      ij = ij + 1
      i1 = ij
      do 11 j = i, n
      zz = h(ij)
      h(ij) = h(ij) / z
      jk = ij
      ik = i1
      do 12 k = i, j
      jk = jk + np - k
      h(jk) = h(jk) - h(ik) * zz
      ik = ik + 1
  12  continue
      ij = ij + 1
  11  continue
      if (h(ij) .le. 0.d0) return
  15  continue
      ij = np
      dmin = h(1)
      do 16 i = 2, n
      if (h(ij) .ge. dmin) go to 16
      dmin = h(ij)
  16  ij = ij + np - i
      if (dmin .le. 0.d0) return
      z = f
      itn = 0
      call funct (n, x, f,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om, normOrNot)
      ifn = 1
      df = dfn
      if (dfn .eq. 0.d0) df = f - z
      if (dfn .lt. 0.d0) df = abs (df * f)
      if (df .le. 0.d0) df = 1.d0
  17  continue
      do 19 i = 1, n
      w(i) = x(i)
  19  continue
      link = 1
      if (idiff - 1) 100, 100, 110
  18  continue
      if (ifn .ge. maxfn) go to 90
  20  continue
  21  continue
      itn = itn + 1
      w(1) = -g(1)
      do 22 i = 2, n
      ij = i
      i1 = i - 1
      z = -g(i)
      do 23 j = 1, i1
      z = z - h(ij) * w(j)
      ij = ij + n - j
  23  continue
  22  w(i) = z
      w(is+n) = w(n) / h(nn)
      ij = nn
      do 25 i = 1, n1
      ij = ij - 1
      z = 0.d0
      do 26 j = 1, i
      z = z + h(ij) * w(is+np-j)
      ij = ij - 1
  26  continue
  25  w(is+n-i) = w(n-i) / h(ij) - z
      z = 0.d0
      gs0 = 0.d0
      do 29 i = 1, n
      if (z * xm(i) .ge. abs (w(is+i))) go to 28
      z = abs (w(is+i)) / xm(i)
  28  gs0 = gs0 + g(i) * w(is+i)
  29  continue
      aeps = eps / z
      iexit = 2
      if (gs0 .ge. 0.d0) go to 92
      alpha = -2.d0 * df / gs0
      if (alpha .gt. 1.d0) alpha = 1.d0
      ff = f
      tot = 0.d0
      int = 0
      iexit = 1
  30  continue
      if (ifn .ge. maxfn) go to 90
      do 31 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  31  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om, normOrNot)
      ifn = ifn + 1
      if (f1 .ge. f) go to 40
      f2 = f
      tot = tot + alpha
  32  continue
      do 33 i = 1, n
      x(i) = w(i)
  33  continue
      f = f1
      if (int - 1) 35, 49, 50
  35  continue
      if (ifn .ge. maxfn) go to 90
      do 34 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  34  continue
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax,
     $     G0w,G0wand,om, normOrNot)
      ifn = ifn + 1
      if (f1 .ge. f) go to 50
      if ((f1 + f2 .ge. f + f) .and.
     $  (7.0d0 * f1 + 5.0d0 * f2 .gt. 12.0d0 * f)) int = 2
      tot = tot + alpha
      alpha = 2.d0 * alpha
      go to 32
  40  continue
      if (alpha .lt. aeps) go to 92
      if (ifn .ge. maxfn) go to 90
      alpha = half * alpha
      do 41 i = 1, n
      w(i) = x(i) + alpha * w(is+i)
  41  continue
      call funct (n, w, f2,epsk,tpar,ns,piob,Xi,Iwmax, 
     $     G0w,G0wand,om, normOrNot)
      ifn = ifn + 1
      if (f2 .ge. f) go to 45
      tot = tot + alpha
      f = f2
      do 42 i = 1, n
      x(i) = w(i)
  42  continue
      go to 49
  45  continue
      z = 0.1d0
      if (f1 + f .gt. f2 + f2)
     $  z = 1.d0 + half * (f - f1) / (f + f1 - f2 - f2)
      if (z .lt. 0.1d0) z = 0.1d0
      alpha = z * alpha
      int = 1
      go to 30
  49  continue
      if (tot .lt. aeps) go to 92
  50  continue
      alpha = tot
      do 56 i = 1, n
      w(i) = x(i)
      w(ib+i) = g(i)
  56  continue
      link = 2
      if (idiff - 1) 100, 100, 110
  54  continue
      if (ifn .ge. maxfn) go to 90
      gys = 0.d0
      do 55 i = 1, n
      w(i) = w(ib+i)
      gys = gys + g(i) * w(is+i)
  55  continue
      df = ff - f
      dgs = gys - gs0
      if (dgs .le. 0.d0) go to 20
      link = 1
      if (dgs + alpha * gs0 .gt. 0.d0) go to 52
      do 51 i = 1, n
      w(iu + i) = g(i) - w(i)
  51  continue
      sig = 1.d0 / (alpha * dgs)
      go to 70
  52  continue
      zz = alpha / (dgs - alpha * gs0)
      z = dgs * zz - 1.d0
      do 53 i = 1, n
      w(iu+i) = z * w(i) + g(i)
  53  continue
      sig = 1.d0 / (zz * dgs * dgs)
      go to 70
  60  continue
      link = 2
      do 61 i = 1, n
      w(iu+i) = w(i)
  61  continue
      if (dgs + alpha * gs0 .gt. 0.d0) go to 62
      sig = 1.d0 / gs0
      go to 70
  62  continue
      sig = -zz
  70  continue
      w(iv+1) = w(iu+1)
      do 71 i = 2, n
      ij = i
      i1 = i - 1
      z = w(iu+i)
      do 72 j = 1, i1
      z = z - h(ij) * w(iv+j)
      ij = ij + n - j
  72  continue
      w(iv+i) = z
  71  continue
      ij = 1
      do 75 i = 1, n
      z = h(ij) + sig * w(iv+i) * w(iv+i)
      if (z .le. 0.d0) z = dmin
      if (z .lt. dmin) dmin = z
      h(ij) = z
      w(ib+i) = w(iv+i) * sig / z
      sig = sig - w(ib+i) * w(ib+i) * z
      ij = ij + np - i
  75  continue
      ij = 1
      do 80 i = 1, n1
      ij = ij + 1
      i1 = i + 1
      do 80 j = i1, n
      w(iu+j) = w(iu+j) - h(ij) * w(iv+i)
      h(ij) = h(ij) + w(ib+i) * w(iu+j)
      ij = ij + 1
  80  continue
      go to (60, 20), link
  90  continue
      iexit = 3
      go to 94
  92  continue
      if (idiff .eq. 2) go to 94
      idiff = 2
      go to 17
  94  continue
      return
 100  continue
      do 101 i = 1, n
         z = hh * xm(i)
         w(i) = w(i) + z
         call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax, 
     $        G0w,G0wand,om, normOrNot)
      g(i) = (f1 - f) / z
      w(i) = w(i) - z
 101  continue
      ifn = ifn + n
      go to (18, 54), link
 110  continue
      do 111 i = 1, n
      z = hh * xm(i)
      w(i) = w(i) + z
      call funct (n, w, f1,epsk,tpar,ns,piob,Xi,Iwmax , 
     $     G0w,G0wand,om, normOrNot)
      w(i) = w(i) - z - z
      call funct (n, w, f2, epsk,tpar,ns,piob,Xi,Iwmax, 
     $     G0w,G0wand,om, normOrNot)
      g(i) = (f1 - f2) / (2.d0 * z)
      w(i) = w(i) + z
 111  continue
      ifn = ifn + n + n
      go to (18, 54), link
      end 
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine difference(nbparm,x,f,epsk,tpar,ns,piob,Xi,
     $     Iwmax,G0w,G0wand,om, normOrNot)
        implicit none
        integer ns,Iwmax,nbparm,icount,i,nsym
        real*8  tpar(ns),epsk(ns),piob,diff,f
        complex*16  cdummy1,G0w(0:Iwmax),G0wand(0:Iwmax)
        complex*16 Xi
        real*8  x(nbparm),om(0:Iwmax),norm
        logical symm, normOrNot
        common /symm/ symm

        icount=1

        if(symm) then
           nsym=(ns-1)/2
           do i=1,nsym
              icount=icount+1
              Epsk(icount)=x(i)
           end do
           icount=0
           do i=nsym+1,2*nsym
              icount=icount+1
              tpar(icount)=x(i)
           end do
        
           do i=nsym+2,ns
              Epsk(i)=-Epsk(i-nsym)
              tpar(i-1)=tpar(i-nsym-1)
           enddo 

            if (normOrNot) then
               norm=0.d0
               
               do i=1,ns-1
                  norm=norm+tpar(i)**2                                            !sum_l V_l^2 = norm
               enddo
               norm=dsqrt(0.25d0/norm)

               do i=1,ns-1
                  tpar(i)=norm*tpar(i)
               enddo
            endif
           

           
        else
      
           do i=1,ns-1
              icount=icount+1
              Epsk(icount)=x(i)
           end do
           icount=0
           do i=ns,2*(ns-1)
              icount=icount+1
              tpar(icount)=x(i)
           end do

           if (normOrNot) then
               norm=0.d0
               
               do i=1,ns-1
                  norm=norm+tpar(i)**2                                            !sum_l V_l^2 = norm
               enddo
               norm=dsqrt(0.25d0/norm)

               do i=1,ns-1
                  tpar(i)=norm*tpar(i)
               enddo
            endif

        endif


        diff=0.d0
        do i=0,Iwmax
           call calcg0(om(i),cdummy1,tpar,epsk,ns,Xi,Iwmax)
           g0wand(i)=cdummy1
c           diff = diff + (abs(g0w(i)-g0wand(i))/om(i))**2
           diff = diff + abs(g0w(i)-g0wand(i))/dfloat(i+1)
c           diff = diff + abs(g0w(i)-g0wand(i))!/dfloat(i+1)
        end do
        f=diff/dfloat(Iwmax+1)

       
        if(symm) then
          
           icount=1
           do i=1,nsym
              icount=icount+1
              x(i)=Epsk(icount)
           end do
           icount=0
           do i=nsym+1,2*nsym
              icount=icount+1
              x(i)=tpar(icount)
           end do
        else
           icount=1
           do i=1,ns-1
              icount=icount+1
              x(i)=Epsk(icount)
           enddo
           icount=0
           do i=ns,2*(ns-1)
              icount=icount+1
              x(i)=tpar(icount)
           enddo
        endif

        end
C========+=========+=========+=========+=========+=========+=========+=$
        subroutine calcg0(omega,g0and,tpar,epsk,ns,Xi,Iwmax)
        implicit none
        integer ns,Iwmax,i
        real*8 tpar(ns),epsk(ns)
        complex*16 g0and,Xi
        real*8 omega
cc
cc      use simple formula for the G_0 function
cc
        g0and=Xi*omega-Epsk(1)
        do i=1,ns-1
           g0and=g0and-tpar(i)**2/(Xi*omega-Epsk(1+i))
        end do

        g0and=dconjg(g0and)/(dconjg(g0and)*g0and)
        end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine datain(nss,uhub,xmu,hmag,ns,beta,wlow,wup,
     $     deltino,imaxmu,deltamu,iteramax, testdmft,dens_wanted,ifix,
     $     inew,iexp,th0,iauto)
c
      implicit none
c     hopping 't', Hubbard U
      real*8 uhub,xmu,hmag,deltamu
      real*8 wlow,wup,deltino,aaa
      real*8 beta,testdmft,dens_wanted,th0
      integer ns,imaxmu,iteramax,nss,is,ival,i,ifix,inew,iexp
      integer iauto
      character*80 text,text1
c
c     read(30,'(a)') text
c     read(30,'(a)') text1
      read(30,'(a)') text
c     hamiltonian parameters
      read(30,*) uhub,hmag
      read(30,'(a)') text
      read(30,*) beta, wlow, wup, deltino
c      read(45,*) varbeta
       
     

c*    
c*==> print options
c*
c*==> max. no. of lines of hamiltonian to print
c*    

      read(30,'(a)') text
      read(30,*) ns,imaxmu,deltamu,iteramax,testdmft   
      read(30,'(a)') text
      read(30,*)ifix,dens_wanted,inew,iauto
      read(30,'(a)') text
      read(30,*)th0,iexp
      write(6,69)ns
 69   format(1x,'# of sites      = ',i5,/)
      if (nss.lt.ns) then
         print*,'nss has to be at least  ',ns
         stop
      endif
      write(6,*)'number of conduction bath levels =',ns-1
      write(6,*)'number of iterations :',iteramax

      call rheader(15)
      read(15,*)
      
      do i=2,ns
         read(15,*)
      end do
      
      read(15,*)
      
      do i=1,ns-1
         read(15,*)
      end do
      
      read(15,*) xmu 


      return
      end
c*    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function isval(b,ns,abas,nss)
c*
      integer b(nss)
      integer *1 abas(nss)
c*
      isval = 0
      do i=1,ns
         isval = isval + (abas(i)+1) * b(i)
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      integer function nchoos(n1,n2)
      xh = 1.0
      if(n2.lt.0) then
        nchoos = 0
        return
      endif
      if(n2.eq.0) then
        nchoos = 1
        return
      endif
      do 100 i = 1, n2
        xh = xh * float(n1+1-i)/float(i)
 100  continue
      nchoos = int(xh + 0.5)
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine statel(ival,b,nss,ns,anew)
      implicit none
      integer nss,ns,ival,i,is,nb
      integer *1 anew(nss)
      integer b(nss)
      is = ival
      do 1 i = ns,1,-1
        nb = is/b(i)
        anew(i) = nb - 1
        is = is - nb * b(i)
    1 continue
c     print*,'anew =',anew
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine computematrix(ioffset,nblock,ntot,b,ns,abas,
     $     nsp1,nmax,nleng,imat,idmat,iddomat,ntable,idouble)
c      program hamiltoniana
c     this is meant to be a subroutine that generates
c     1) the blocks
c     2) the non-diagonal matrix elements within each block (imat)
c     3) the matrix elements of d^{dagger} -> idmat
      implicit none

      integer nblock,i,nup,ndo,nchoos,j,nmax
      integer ns,nsp1,ntot,info
      integer ioffset(nblock+1)
      integer b(ns)
      integer ntable(ntot)
      integer*1 abas(ns)
      integer nleng(nblock)
      integer*1 imat(nmax,nmax,nblock)
      integer idmat(ntot)
      integer iddomat(ntot)
      integer*1 idouble(ntot)
      
      
      
c     compute ioffsets
      
      ioffset(1)=0
      ioffset(2)=1
      do i=2,nblock-1
         nup=mod(i-1,nsp1)
         ndo=(i-1)/nsp1
         ioffset(i+1)=ioffset(i)+nchoos(ns,nup)*nchoos(ns,ndo)
c         write(6,*) 'i', i, 'up', nup, 'down', ndo, 'ioffsets',ioffset(i)
      enddo
      ioffset(nblock+1)=ntot
      
      call buildblocks(ns,ntot,nblock,ntable,ioffset,nleng,
     $     b,abas,nsp1)
      
c     loop over blocks -> compute the Hamiltonian

      do i=1,ntot
         idouble(i)=0
      enddo

      do i=1,nblock
c         write(6,*) i, nleng(i)
         do j=1,nleng(i)
            call findnupndo(ntable(ioffset(i)+j)-1,nup,ndo,b,
     $           abas,ns)
            if (abas(1).eq.2) idouble(ioffset(i)+j)=1
         enddo
         call hamilt(ntable(ioffset(i)+1),nleng(i),
     $        imat(1,1,i),b,ns,abas,nmax)
         if (mod(i-1,nsp1).ne.ns) then
            call ddag(ntable(ioffset(i)+1),ntable(ioffset(i+1)+1),
     $           nleng(i),nleng(i+1),idmat(ioffset(i)+1),b,ns,
     $           abas,ntot,i)
         endif
         if(i.gt.(ns+1)) then
            call ddown(ntable(ioffset(i)+1),ntable(ioffset(i-ns-1)+1),
     $           nleng(i),nleng(i-ns-1),iddomat(ioffset(i)+1),b,ns,
     $           abas,ntot,i)
         endif
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine hamilt(ntablein,nlen,imatin,b,ns,abas,nmax)
      implicit none
      integer nlen,ns,nmax,i,j,nup,ndo,ntnew,inew,k,isegno
      integer ntablein(nlen)
      integer*1 imatin(nmax,nmax)
      integer b(ns)
      integer*1 abas(ns)

      do i=1,nlen
         do j=1,nlen
            imatin(i,j)=0
         enddo
      enddo

  
c     diagonal part is already coded in ntab!!!!
c     non-diagonal part

      do i=1,nlen
         call findnupndo(ntablein(i)-1,nup,ndo,b,abas,ns)
c     if up is on site 1, it can hop
c         write(6,*) ntablein(i)
         if (abas(1).gt.0) then
c     loop over final sites
            do j=2,ns
               if (abas(j).le.0) then
                  ntnew=ntablein(i)+b(j)*(1-2*abas(j))+b(1)*
     $                 (1-2*abas(1))
                  call findinew(ntnew,ntablein,nlen,inew)
                  nup=0
                  do k=2,j-1
                     if (abas(k).gt.0) nup=nup+1
                  enddo
                  isegno= 1-2*mod(nup,2)
                  imatin(i,inew)=j*isegno
                  imatin(inew,i)=j*isegno
c                  write(6,*) 'inew', inew, 'imatin', imatin(i,inew)
               endif
            enddo
         endif
c     if down is on site 1, it can hop
         if (abas(1).eq.2.or.abas(1).eq.-1) then
            do j=2,ns
               if (abas(j).eq.0.or.abas(j).eq.1) then
                  ntnew=ntablein(i)+b(j)*(2*abas(j)-1)-b(1)*
     $                 abas(1)/abs(abas(1))
                  call findinew(ntnew,ntablein,nlen,inew)
                  
                  ndo=0
                  do k=2,j-1
                     if (abas(k).eq.2.or.abas(k).eq.-1) ndo=ndo+1
                  enddo
                  isegno=1-2*mod(ndo,2)
                  imatin(i,inew)=j*isegno
                  imatin(inew,i)=j*isegno                  
               endif
            enddo
         endif
        enddo


      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findinew(ntnew,ntab,nlen,inew)
      implicit none
      integer nlen,i,ntnew,inew
      integer ntab(nlen)
      do i=1,nlen
         if (ntab(i).eq.ntnew) goto 777
      enddo

      write(6,*)'il numero ',ntnew,' manca'
      
 777  inew=i
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine buildblocks(ns,ntot,nblock,ntable,ioffset,nfound,
     $     b,abas,nsp1)
      implicit none

      integer ns,ntot,nblock
      integer ntable(ntot)
      integer ioffset(nblock+1)
      integer nfound(nblock)
      integer b(ns)
      integer*1 abas(ns)
      integer nup,ndo,nsp1,indblock,i

      do i=1,ntot
         call findnupndo(i-1,nup,ndo,b,abas,ns)
         indblock=nup+nsp1*ndo+1
         nfound(indblock)=nfound(indblock)+1

         if (nfound(indblock).gt.(ioffset(indblock+1)-
     $        ioffset(indblock)))
     $        then
            write(6,*)'Il blocco ',indblock,
     $           ' deve essere lungo', nfound(indblock)
            write(6,*)'nup, ndo =',nup,ndo
            write(6,*)'ci voglio mettere',i
            stop
         endif
         ntable(nfound(indblock)+ioffset(indblock))=i
c         write(6,*) ntable
      enddo
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine findnupndo(nstato,nup,ndo,b,anew,ns)
      implicit none


      integer nb,i
      integer ns,nup,ndo,nstato,is
      integer *1 anew(ns)
      integer b(ns)

      is = nstato
      nup=0
      ndo=0
      do 1 i = ns,1,-1
        nb = is/b(i)
        anew(i) = nb - 1
        is = is - nb * b(i)
        if (anew(i).gt.0) nup=nup+1
        if (anew(i).eq.2.or.anew(i).eq.-1) ndo=ndo+1
    1 continue
      write(72,*)(nstato+1),(anew(i),i=1,ns)
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ddag(ntab0,ntab1,
     $     nlen0,nlen1,idmatin,b,ns,abas,ntot,iblocco)
      implicit none

      integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k
      integer ntab0(nlen0)
      integer ntab1(nlen1)
      integer b(ns)
      integer*1 abas(ns)
      integer ntot,iblocco
      integer idmatin(nlen0)
      do j=1,nlen0
         iold=ntab0(j)
         call findnupndo(iold-1,nup,ndo,b,abas,ns)
         if (abas(1).le.0) then
            inew=iold+(1-2*abas(1))
            do k=1,nlen1
               if (ntab1(k).eq.inew) then
                  idmatin(j)=k
                  goto 555
               endif
            enddo
            write(6,*)'Non trovo lo stato ottenuto applicando d+'
 555        continue
         else
            idmatin(j)=0
         endif
      enddo
      return
      end

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ddown(ntab0,ntab1,
     $     nlen0,nlen1,iddomatin,b,ns,abas,ntot,iblocco)
      implicit none

      integer nlen0,nlen1,ns,j,iold,nup,ndo,inew,k,isegno
      integer ntab0(nlen0)
      integer ntab1(nlen1)
      integer b(ns)
      integer*1 abas(ns)
      integer ntot,iblocco
      integer iddomatin(nlen0)
      do j=1,nlen0
         iold=ntab0(j)
         call findnupndo(iold-1,nup,ndo,b,abas,ns)


cccc     If down is on site 1

 
         if (abas(1).eq.-1.or.abas(1).eq.2) then
            inew=iold + (1-2*abas(1))/3

ccccccccccc          check delle operazioni   ccccccccccccccccc
 
c            write(6,*) 'abas', abas
c            write(6,*)'in', inew,'io', iold,'dsum',(1-2*abas(1))/3
c            call findnupndo(inew-1,nup,ndo,b,abas,ns)
c            write(6,*) 'abasnew', abas
c            call findnupndo(inew-1,nup,ndo,b,abas,ns)
c            write(6,*) 'inew', inew, 'abas', abas
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  

cc          ddown operator has to hop nup spin to reach his spin down
          
            isegno=1-2*mod(nup,2)
       
            do k=1,nlen1
               if (ntab1(k).eq.inew) then
                  iddomatin(j)=k*isegno
                  goto 555
               endif
            enddo
            

            write(6,*)'Non trovo lo stato ottenuto applicando d'
 555        continue
         else
            iddomatin(j)=0
         endif
      enddo
      return
      end


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diag(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,          !Integer-Input   
     $                blocksize,
     $                b,ntable,ioffset,imat,                           !Integer-Array-Input
     $                xmu,uhub,hmag,                                   !Real-Input
     $                tpar,epsk,                                       !Real-Array-Input
     $                emin,                                            !Real-Output
     $                eig,zeig)                                        !Real-Array-Output

      implicit none

!Input-Variables
      integer ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,blocksize
      integer b(ns),ntable(ntot),ioffset(nblock+1)
      integer*1 imat(nmaxx,nmaxx)
      real*8 xmu,uhub,hmag
      real*8 tpar(ns),epsk(ns)
!Output-Variables
      real*8 emin
      real*8 eig(nmaxx),zeig(nmaxx,nmaxx)
!Subroutine-Internal-Variables
      integer j,k,is,neig,nup,ndo,icount,info
      integer*1 abas(ns)
      real*8 work(3*nmaxx)
      real*8 realmat(nmaxx*(nmaxx+1)/2),xmat(nmaxx,nmaxx)
      
      epsk(1)=-xmu
      emin=1.d16

      do j=1,nmaxx
         do k=j+1,nmaxx
            xmat(j,k)=0.d0
            xmat(k,j)=0.d0
         enddo
         xmat(j,j)=100000.d0
      enddo

      do j=1,blocksize
         call findnupndo(ntable(ioffset(myblock)+j)-1,nup,ndo,
     $        b,abas,ns)
         xmat(j,j)=0.d0
         do is=1,ns
            xmat(j,j)=xmat(j,j)+epsk(is)*dabs(dfloat(abas(is)))
         enddo
         if (abas(1).eq.2) then
            xmat(j,j)=xmat(j,j)+uhub
         else
            xmat(j,j)=xmat(j,j)-hmag*dfloat(abas(1))
         endif
         do k=j+1,blocksize
            if (imat(j,k).ne.0) then
               xmat(j,k)=xmat(j,k)+tpar(abs(imat(j,k))-1)*
     $                             dfloat(imat(j,k))/
     $                             abs(dfloat(imat(j,k)))
               xmat(k,j)=xmat(k,j)+tpar(abs(imat(j,k))-1)*
     $                             dfloat(imat(j,k))/
     $                             abs(dfloat(imat(j,k)))
            endif
         enddo
      enddo
         
      icount=0
         
      do j=1,nmaxx
         do k=j,nmaxx
            icount=icount+1
            realmat(icount)=xmat(k,j)
         enddo
      enddo
            
       
      call dspev('V','L',nmaxx,realmat,eig,zeig,
     $           nmaxx,work,info)

      do neig=1,blocksize
         if (eig(neig).lt.emin) emin=eig(neig)
      enddo

      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine lehmann(ns,nsp1,nblock,ntot,nmax,nmaxx,myblock, !Integer-Input   
     $                   blocksize,blocksize1,
     $                   ioffset,idouble,idmat,                  !Integer-Array-Input
     $                   xmu,uhub,hmag,beta,                     !Real-Input
     $                   eig,zeig,eigold,zold,                   !Real-Array-Input
     $                   double,zpart,dens,                      !Real-Output
     $                   rlehm)                                  !Real-Array-Output

      implicit none
!Input-Variables
      integer ns,nsp1,nblock,ntot,nmax,nmaxx,myblock,blocksize
      integer blocksize1
      integer ioffset(nblock+1),idmat(ntot)
      integer*1 idouble(ntot)
      real*8 xmu,uhub,hmag,beta
      real*8 eig(nmaxx),zeig(nmaxx,nmaxx)
      real*8 eigold(nmaxx),zold(nmaxx,nmaxx)
!Output-Variables
      real*8 double,zpart,dens
      real*8 rlehm(nmaxx,nmaxx)
!Subroutine-Internal-Variables
      integer i,j,k,meig,neig
      real*8 rmel

      double=0.0d0
      zpart=0.0d0
      dens=0.0d0

      do i=1,nmaxx
         do j=1,nmaxx
            rlehm(i,j)=0.d0
         enddo
      enddo

      if (mod(myblock-1,nsp1).ne.0) then
c     loop over |ndo,nup-1> basis states
         do neig=1,blocksize
            do meig=1,blocksize1
               rmel=0.d0
               do j=1,blocksize1
                  if (idmat(ioffset(myblock-1)+j).ne.0) then
                     rmel=rmel+zold(j,meig)*
     $                         zeig(idmat(ioffset(myblock-1)+j),neig)
                  endif
               enddo
               rlehm(neig,meig)=rmel**2
               dens=dens+rlehm(neig,meig)*dexp(-beta*eig(neig))
            enddo
         enddo
      endif
      
      dens=2.0d0*dens

      do neig=1,blocksize
         do j=1,blocksize
            double=double+
     $           zeig(j,neig)*zeig(j,neig)*
     $           dfloat(idouble(ioffset(myblock)+j))*
     $           dexp(-beta*eig(neig))
         enddo
         zpart=zpart+dexp(-beta*eig(neig))
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine computegimag(omoffset,omnumber,                 !Integer-Input
     $                        nsp1,nmaxx,Iwmax,prozessoren,         
     $                        nleng,                             !Integer-Array-Input  
     $                        beta,                              !Real-Input
     $                        eig,rlehm,om,                      !Real-Array-Input
     $                        Xi,                                !Complex-Input
     $                        Gw)                                !Complex-Array-Output

      implicit none

!Input-Variables
      integer omoffset,omnumber,nsp1,nmaxx,Iwmax,prozessoren
      integer nleng(nsp1**2)
      real*8 beta
      real*8 eig(nmaxx,prozessoren)
      real*8 rlehm(nmaxx,nmaxx,prozessoren),om(0:Iwmax)
      complex*16 Xi
!Output-Variables
      complex*16 Gw(0:Iwmax)
!Subroutine-Internal-Variables
      integer i,iomega,meig,neig
      
      do i=0,Iwmax
         Gw(i)=dcmplx(0.d0,0.d0)
      enddo

      do i=1,nsp1**2
      if (mod(i-1,nsp1).ne.0) then
         do iomega=omoffset+1,omoffset+omnumber
            do neig=1,nleng(i)
               do meig=1,nleng(i-1)
                  Gw(iomega)=Gw(iomega)+rlehm(neig,meig,i)
     $                 /(Xi*om(iomega)-
     $                 (eig(neig,i)-eig(meig,i-1)))*
     $                 (dexp(-beta*eig(neig,i))+
     $                  dexp(-beta*eig(meig,i-1)))
               enddo                  
            enddo
         enddo
      endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine computegreal(omoffset,omnumber,                 !Integer-Input
     $                        nsp1,nmaxx,Iwmaxreal,prozessoren,         
     $                        nleng,                             !Integer-Array-Input  
     $                        beta,                              !Real-Input
     $                        eig,rlehm,                         !Real-Array-Input
     $                        omr,                               !Complex-Arral-Input
     $                        Gw)                                !Complex-Array-Output

      implicit none

!Input-Variables
      integer omoffset,omnumber,nsp1,nmaxx,Iwmaxreal,prozessoren
      integer nleng(nsp1**2)
      real*8 beta
      real*8 eig(nmaxx,prozessoren)
      real*8 rlehm(nmaxx,nmaxx,prozessoren)
      complex*16 omr(0:Iwmaxreal)
!Output-Variables
      complex*16 Gw(0:Iwmaxreal)
!Subroutine-Internal-Variables
      integer i,iomega,meig,neig
      
      do i=0,Iwmaxreal
         Gw(i)=dcmplx(0.d0,0.d0)
      enddo
      
      do i=1,nsp1**2
      if (mod(i-1,nsp1).ne.0) then
         do iomega=omoffset+1,omoffset+omnumber
            do neig=1,nleng(i)
               do meig=1,nleng(i-1)
                  Gw(iomega)=Gw(iomega)+rlehm(neig,meig,i)
     $                 /(omr(iomega)-
     $                 (eig(neig,i)-eig(meig,i-1)))*
     $                 (dexp(-beta*eig(neig,i))+
     $                  dexp(-beta*eig(meig,i-1)))
               enddo                  
            enddo
         enddo
      endif
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine dos2(omnumber,omoffset,Iwmaxreal,ksteps,        !Integer-Input
     $                beta,xmu,Pi,                               !Real-Input  
     $                Energya,                                   !Real-Array-Input
     $                omr,sigre,                                 !Complex-Array-Input
     $                dos)                                       !Real-Array-Output

      implicit none
     
!Input-Variables
      integer omnumber,omoffset,Iwmaxreal,ksteps
      real*8 beta,xmu,Pi
      real*8 Energya(0:ksteps,0:ksteps)
      complex*16 omr(0:Iwmaxreal),sigre(0:Iwmaxreal)
!Output-Variables
      real*8 dos(0:Iwmaxreal)
!Subroutine-Internal-Variables
      integer i,ikx,iky
      real*8 weightx,weighty,weightm
 
      if (omoffset.eq.-1)then
         write(6,*) 'Computing the interacting DOS'
         call flush(6)
      endif

      do i=omoffset+1,omoffset+omnumber
         dos(i)=0.d0
         do ikx=0,ksteps
            weightx=1.0d0
            weightm=1.0d0
            if ((ikx.eq.0).or.(ikx.eq.ksteps)) then
              weightx=0.5d0
            endif
            do iky=0,ikx
               weighty=1.0d0
               if ((iky.eq.0).or.(iky.eq.ikx)) then
                  weighty=0.5d0
               endif
               if (((iky.eq.0).and.(ikx.eq.0)).or.
     $             ((iky.eq.ksteps).and.(ikx.eq.ksteps))) then
                  weightm=0.5d0
               endif
               dos(i)=dos(i)-(2.0d0*weightx*weighty*weightm/Pi)*
     $                dimag(1.d0/(omr(i)-Energya(ikx,iky)+xmu
     $                -sigre(i)))
            enddo
         enddo
         dos(i)=dos(i)/(dfloat(ksteps)*dfloat(ksteps))
      enddo

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine selfconst2(omnumber,omoffset,Iwmax,ksteps,      !Integer-Input
     $                      beta,xmu,                            !Real-Input  
     $                      om,Energya,                          !Real-Array-Input
     $                      Xi,                                  !Complex-Input      
     $                      g0wand,Gw,                           !Complex-Array-Input
     $                      Gloc)                                !Complex-Array-Output
      
      implicit none

!Input-Variables
      integer omnumber,omoffset,Iwmax,ksteps,myblock
      real*8 beta,xmu
      real*8 om(0:Iwmax),Energya(0:ksteps,0:ksteps)
      complex*16 Xi
      complex*16 g0wand(0:Iwmax),Gw(0:Iwmax)
!Output-Variables
      complex*16 Gloc(0:Iwmax)
!Subroutine-Internal-Variables
      integer i,ikx,iky
      real*8 weightx,weighty,weightm
      complex*16 W(omoffset+1:omoffset+omnumber)
      complex*16 gand(omoffset+1:omoffset+omnumber)

      do i=omoffset+1,omoffset+omnumber
         gand(i)=1.d0/g0wand(i)
         Gw(i)=1.d0/Gw(i)
         W(i)= Xi*om(i)-gand(i)+Gw(i)+xmu!-tpri
      enddo 
      
      do i=omoffset+1,omoffset+omnumber
         Gloc(i)=dcmplx(0.d0,0.d0)
         do ikx=0,ksteps
            weightx=1.0d0
            weightm=1.0d0
            if ((ikx.eq.0).or.(ikx.eq.ksteps)) then
              weightx=0.5d0
            endif
            do iky=0,ikx
              weighty=1.0d0
              if ((iky.eq.0).or.(iky.eq.ikx)) then
                weighty=0.5d0
              endif
              if (((iky.eq.0).and.(ikx.eq.0)).or.
     $            ((iky.eq.ksteps).and.(ikx.eq.ksteps))) then
                 weightm=0.5d0
              endif
              Gloc(i)=Gloc(i)+2.0d0*weightx*weighty*weightm/(W(i)-
     $             Energya(ikx,iky))
            enddo
         enddo
         Gloc(i)=Gloc(i)/(dfloat(ksteps)*dfloat(ksteps))
         Gloc(i)=1.d0/Gloc(i)+gand(i)-Gw(i)
         Gloc(i)=0.5d0/Gloc(i)+0.5d0*g0wand(i)
      enddo
   
      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine selfconst2_magnetic(omnumber,omoffset,Iwmax,ksteps,L,    !Integer-Input
     $                               beta,xmu,B,                          !Real-Input  
     $                               om,Energya,                          !Real-Array-Input
     $                               Xi,                                  !Complex-Input      
     $                               g0wand,Gw,                           !Complex-Array-Input
     $                               gloc)                                !Complex-Array-Output
      
      implicit none

!Input-Variables
      integer omnumber,omoffset                                             !omnumber is the number of all omegas, omoffset normally -omnumber/2
      integer Iwmax                                                         !maximum value of omega frequency
      integer ksteps                                                        !twice many steps in the kGrid, e.g. ksteps = 2 gives kx = -pi/2, 0, pi/2, pi
      integer kstepsX
      integer kstepsY                                                       !since Brillouin in y direction smaller but same finesse, have less ksteps
      integer myblock
      integer L                                                             !denominator in B, in Georgs paper it is q
      real*8 beta                                                           !Beta as 1/T with temperature T
      real*8 xmu                                                            !chemical potential
      real*8 B                                                              !magnetic field is given from outside and can have p as well
!B...magnetic field in units of [Phi_0 (Fluxquantum)/area of unit cell(=a^2=1)]
!L...Size of unit cell in y-direction -> normally choosen such that
!for testing: set L finite and B=0!
      real*8 om(0:Iwmax),Energya(0:ksteps,0:ksteps)                         !om are the frequencies omega, Energya are tight binding energies
      complex*16 Xi                                                         !just the imaginary number i
      complex*16 g0wand(0:Iwmax),Gw(0:Iwmax)                                !g0wand is G_0 of anderson Gw is full Greens function
!Output-Variables
      complex*16 gloc(0:Iwmax,1:L,1:L)                                      !the output function which is the local Greens function G_{L, L^'}(omega) of each band L in dependence of omega 
!Subroutine-Internal-Variables
      integer i,j,k,ikx,iky                                                 !running variables in do loops
      real*8 weightx,weighty,weightm                                        !since kx and ky are symmetric to 0 we do not have to calculate all, but can multiply them with weights, eg 4 at the diagonale
      real*8 dkx, dky
      real*8 kx,ky                                                          !real values of kx and ky in a Brillouin zone and not the running value ikx and iky
      real*8 Pi                                                             !will store the value of pi
      real*8 t,t1,t2,checkone                                               !t is hopping parameter, checkone to check the momentum sum

      complex*16 W(omoffset+1:omoffset+omnumber)                            !basically i\nu + \mu - (G_0^{-1} - G^{-1}) = i\nu + \mu - \Sigma
      complex*16 gand(omoffset+1:omoffset+omnumber)                         !anderson Greens function
      complex*16 Gw1(omoffset+1:omoffset+omnumber)                          !G^{-1}
      complex*16 epsmat(1:L,1:L)                                            !the epsilon matrix I give
      complex*16 ginv(1:L,1:L)                                              !only used to store G^{-1}
      complex*16 gk                                                         !never needed
!Varaibles for inversion subroutines
      INTEGER infoinv                                                       !written by Lapack, if exited succesfull or not
      INTEGER ipiv(1:L)                                                     !contaions some pivot indices for Lapack
      COMPLEX*16 workinv(1:10*L)                                            !only needed in Lapack and has those dimensions, can return optimal WORK
      include 'tpri.dat'
      Pi=dacos(-1.0d0)                                                      !get value of Pi by the arccos(-1)

!      write(6,*)B

      do i=omoffset+1,omoffset+omnumber                                    !We basically go from -omnumber + 1 to omnumber and therefore we look at each frequency by using the running variable of omega
         gand(i)=1.d0/g0wand(i)                                             !G_0^{-1}
         Gw1(i)=1.d0/Gw(i)                                                  !G^{-1}
         W(i)= Xi*om(i)-gand(i)+Gw1(i)+xmu!-tpri                            !basically i\nu + \mu - (G_0^{-1} - G^{-1}) = i\nu + \mu - \Sigma
      enddo

!Initialize dispersion matrix
      do j=1,L
         do k=1,L
            epsmat(j,k)=dcmplx(0.0d0,0.0d0)                             !Set each entry of the epsilon matrix 0
            do i=0,Iwmax
               gloc(i,j,k)=dcmplx(0.0d0,0.0d0)                          !Set each entry in G_{L, L^'}(omega) 0
            enddo
         enddo
      enddo

      checkone=0.0d0
      kstepsX = ksteps
      kstepsY = ksteps/L                                                         !Maybe change this .................
      dkx = (2 * Pi)/dfloat(kstepsX)
      dky = (2 * Pi)/dfloat(kstepsY * L)                                 !Make sure to have float division 

      do ikx = 1, kstepsX                                                !get kx from the running variable ikx
        kx = -Pi + dfloat(ikx - 1) * dkx
        do iky = 1, kstepsY
         ky = -Pi/dfloat(L) + dfloat(iky - 1) * dky
         checkone=checkone+1.0d0                                     !by this we sum over all (kx,ky) tuples which should be ksteps * ksteps/L = 2*ksteps * 2*kstepsY using integer division 

!OLD SAMPLING
!do ikx=-ksteps+1,ksteps
!   kx=Pi*dfloat(ikx)/dfloat(ksteps)                               !get kx from the running variable ikx
!   do iky=-kstepsY+1,kstepsY                                      !it is totally fine to have less points in y direction since the finesse should stay the same and ky is shorter in generell
!      ky=Pi*dfloat(iky)/dfloat(ksteps)                            !therefore do not need changes here, since it is the same besides less k points
!      checkone=checkone+1.0d0                                     !by this we sum over all (kx,ky) tuples which should be 2*ksteps * 2*ksteps/L = 2*ksteps * 2*kstepsY using integer division 
   
!Construct dispersion matrix for given values of kx and ky
            if (L.gt.2) then                                            !if epsilon matrix is bigger than a 2x2 matrix, L is bigger than 2
               epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)
               epsmat(1,2)=dcmplx(-t,0.0d0)
               epsmat(1,L)=dcmplx(-t*dcos(dfloat(L)*ky),0.0d0)+         !top right entry
     $              dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
               epsmat(L,L)=dcmplx(-2.0d0*t*                             !bottom right entry
     $              dcos(kx+dfloat(L-1)*B*2.0d0*Pi),0.0d0)
               epsmat(L,L-1)=dcmplx(-t,0.0d0)                           !bottom right but one to the left entry
               epsmat(L,1)=dcmplx(-t*dcos(dfloat(L)*ky),0.0d0)-         !bottom left entry
     $              dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))                 !used cos and sin instead of Euler because maybe faster
               
               do j=2,L-1
                  epsmat(j,j)=dcmplx(-2.0d0*t*                          !diagonale entries
     $                 dcos(kx+dfloat(j-1)*B*2.0d0*Pi),0.0d0)
                  epsmat(j,j-1)=dcmplx(-t,0.0d0)                        !lower diagonale entries
                  epsmat(j,j+1)=dcmplx(-t,0.0d0)                        !upper diagonale entries
               enddo
            else                                                        !if epsilon matrix is a 2x2 matrix
               epsmat(1,1)=dcmplx(-2.0d0*t*dcos(kx),0.0d0)              !top left entry
               epsmat(1,2)=dcmplx(-t-t*dcos(dfloat(L)*ky),0.0d0)+       !top right entry
     $              dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
               epsmat(2,2)=dcmplx(-2.0d0*t*dcos(kx+B*2.0d0*Pi),0.0d0)   !bottom right entry
               epsmat(2,1)=dcmplx(-t-t*dcos(dfloat(L)*ky),0.0d0)-       !bottom left entry
     $              dcmplx(0.0d0,-t*dsin(dfloat(L)*ky))
            endif

            do i=omoffset+1,omoffset+omnumber         
               do j=1,L
                  do k=1,L
                     ginv(j,k)=-epsmat(j,k)                             !do we need i loop here? Still correct but needless but you need less memory
                  enddo
               enddo
               do j=1,L
                  ginv(j,j)=ginv(j,j)+W(i)                              ![-\epsilon(k)_{j,j} + i\nu + \mu - \Sigma] = G_{j,j}^{-1}(\nu, k), but have \epsilon_{l,lPrime} outside diagonale as well
               enddo
            
               CALL ZGETRF(L,L,ginv,L,ipiv,infoinv)                     !Lapack routine for "These subroutines factor general matrix A using Gaussian elimination with partial pivoting"
               CALL ZGETRI(L,ginv,L,ipiv,workinv,10*L,infoinv)          !Lapack routine to "ZGETRI computes the inverse of a matrix using the LU factorization computed by ZGETRF" Now ginv is G(\nu, k), inverts ginv
               
               do j=1,L
                  do k=1,L
                     gloc(i,j,k)=gloc(i,j,k)+ginv(j,k)                  !Now we get Glocal_{j,k}(omega) = sum over k (G_{l,lPrime} (\nu, k))
                  enddo
               enddo
            enddo
         enddo
      enddo 

      !checkone=checkone*dfloat(L)/(4.0d0*dfloat(ksteps**2))
      checkone=checkone/((dfloat(kstepsX))*(dfloat(kstepsY)))                     !if this is not 1.0 or 10^(-15) or so away, then the normalization is wrong
      if (omoffset.le.2) then
         write(6,*)"Check momentum sum: ",checkone                                           !should be exactly 1
      endif

      do i=omoffset+1,omoffset+omnumber                                                      !basically from -omnumber/2 to omnumber/2
         do j=1,L
            do k=1,L
               !gloc(i,j,k)=gloc(i,j,k)*dfloat(L)/(4.0d0*dfloat(ksteps**2))                  !this is the wrong normalization, since it does not account for the integer division when ksteps/L
               gloc(i,j,k)=gloc(i,j,k)/(dfloat(kstepsX)*dfloat(kstepsY))                                   !better do this where we consider integer division by using kstepsY

            enddo
         enddo
      enddo
   
      if (omoffset.le.2) then
         do j=1,L
            do i=omoffset+1,omoffset+omnumber
               write(1000+j,*)om(i),
     $              dreal(gloc(i,j,j)),dimag(gloc(i,j,j))
            enddo
         enddo
      endif


      return 
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine selfconst3(g0wand,Gw,Gloc,Xi,Pi,Iwmax,om,xmu)
    

      integer i,j,k,Iwmax,ix,iy,iz
      real*8 Pi,xmu,Ek,var
      complex*16 Xi
      real*8 om(0:Iwmax)
      complex*16 g0wand(0:Iwmax),Gw(0:Iwmax),Gloc(0:Iwmax)
      complex*16 W(0:Iwmax)
      real*8 Lx,Ly,Lz,Vol,kx,ky,kz,k0x,k0y,k0z,tsc,tfcc,d
      parameter (Lx=100)
      parameter (var=0.5d0)
      parameter(d=4.d0*var)
      

      do i=0,Iwmax
         g0wand(i)=1.d0/g0wand(i)
         Gw(i)=1.d0/Gw(i)
         W(i)= Xi*om(i)+xmu- g0wand(i)+Gw(i)
      enddo 

      
 
      Ly=Lx
      Lz=Lx
      Vol = Lx*Ly*Lz
      k0x = Pi*2.d0/Lx
      k0y = Pi*2.d0/Ly
      k0z = Pi*2.d0/Lz

! fcc/cubic Gitter
! Hopping fcc
c      tfcc = 4.d0/sqrt(6.d0*d*d+12.d0)
! Hopping cubic
c      tsc  = d*tfcc/2.d0
      tsc=var*dsqrt(2.d0)/dsqrt(3.d0)
      write(6,*) 'Variance = ', var, 'D = ',tsc*3.d0
      tfcc=0.d0
! loop over momenta separated into four blocks
! uses symmetry kx -- ky -- kz
! ----------------------------------------------------

      kx = 0.d0
      do ix=1,Lx/2-1
c         if (ix/10*10.eq.ix) write(FOUT,*) 'ix=',ix
         kx = kx+k0x
         ky = 0.d0
         do iy=1,Ly/2-1
            ky = ky+k0y
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
c     here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
c               if (D.le.1.d0) then
c                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
c               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
c               endif
c               rzero=rzero+8.d0
c               rone=rone+8.d0*eq0
c               rtwo=rtwo+8.d0*eq0*eq0
                 
               
c     here comes the loop over Matsubara frequencies
               do i=0,Iwmax
                   Gloc(i)=Gloc(i)+(8.d0,0.d0)/(W(i)-Ek)
                enddo
            enddo
            
         enddo
      enddo
      
! ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = 0.d0
         do iy=1,Ly/2-1
            ky = ky+k0y
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
                                ! here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c               rzero=rzero+12.d0
c               rone=rone+12.d0*eq0
c               rtwo=rtwo+12.d0*eq0*eq0
               
! here comes the loop over Matsubara frequencies
 
               do i=0,Iwmax
                  Gloc(i) =Gloc(i)+(12.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
! ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = -pi
         do iy=1,2
            ky = ky+pi
            kz = 0.d0
            do iz=1,Lz/2-1
               kz = kz+k0z
!     here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek = Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c               rzero=rzero+6.d0
c               rone=rone+6.d0*eq0
c               rtwo=rtwo+6.d0*eq0*eq0
               
! here comes the loop over Matsubara frequencies

               do i=0,Iwmax
                  Gloc(i) = Gloc(i)+ (6.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
c     ----------------------------------------------------

      kx = -pi
      do ix=1,2
         kx = kx+pi
         ky = -pi
         do iy=1,2
            ky = ky+pi
            kz = -pi
            do iz=1,2
               kz = kz+pi
                                ! here comes the relevant dispersion relation
               Ek = tfcc*((cos(kx)+cos(ky))*cos(kz)+cos(kx)*cos(ky))
               if (D.le.1.d0) then
                  Ek=Ek+tsc*(cos(2.d0*kx)+cos(2.d0*ky)+cos(2.d0*kz))
               else  
                  Ek = Ek+tsc*(cos(kx)+cos(ky)+cos(kz))
               endif
c     rzero=rzero+1.d0
c     rone=rone+eq0
c     rtwo=rtwo+eq0*eq0
                 
c     here comes the loop over Matsubara frequencies
      
               do i=0,Iwmax
                  Gloc(i)=Gloc(i)+(1.d0,0.d0)/(W(i)-Ek)
               enddo
               
            enddo
         enddo
      enddo
! ----------------------------------------------------

         


      do i=0,Iwmax
         Gloc(i) = Gloc(i)/Vol
         Gloc(i)=1.d0/Gloc(i)+g0wand(i)-Gw(i)
         Gloc(i)=1.d0/Gloc(i)
      enddo

c         write(FOUT,*)
c         write(FOUT,'(a20,3f16.8)') 'first three moments:'&
c             ,real(rzero/Vol),real(rone/Vol),real(rtwo/Vol)
c         write(FOUT,*)



      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      function EllipK(x1)
c      real*8 x1,a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,EllipK
c      if (x1.le.0.or.x1.ge.1) then
c            print *,x1
c            stop 99
c      endif
c      data a0,a1,a2,a3,a4 /1.38629436112,0.09666344259,0.03590092383,
c     .      0.03742563713,0.01451196212/
c      data b0,b1,b2,b3,b4 /0.5,0.12498593597,0.06880248576,
c     .      0.03328355346,0.00441787012/
c      EllipK=(((a4*x1+a3)*x1+a2)*x1+a1)*x1+a0-
c     .      ((((b4*x1+b3)*x1+b2)*x1+b1)*x1+b0)*dlog(x1)
c      return
c      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      function Energy(x,y)
     
      real*8 x,y
      real*8 Energy,t,t1,t2
      include 'tpri.dat'
 
      Energy=-2.0d0*t*(dcos(x)+dcos(y))-4.0d0*t1*dcos(x)
     $             *dcos(y)-2.0d0*t2*(dcos(2.0d0*x)+dcos(2.0d0*y))
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
