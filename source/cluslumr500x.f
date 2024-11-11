      parameter (ang=3.14159265/180.)
      character fn1*18,fn2*18,fn3*18,aa*100
      real az1(18),az2(10),z1(8000000,7)
      real a5(12),z3(200),z4(200)
      integer id1,id2,ip1,ip2(2),index1(8000000)
      open(12,file='cluslumx.lst')
      open(25,file='bc2003_magnitude_F119')
      N0=0
      do i=1,50
         read(25,'(a)')aa
      enddo
 26    read(25,*,end=27)a5
      N0=N0+1
      z3(N0)=a5(1)
      z4(N0)=a5(12)-28.0
      goto 26
 27   close(25)
      write(*,*)N0
      
 15   read(12,*,end=16)fn1,fn2,fn3
      open(10,file='../data/'//fn1)
      open(20,file=fn2)
      open(30,file=fn3)
      N=0
 11   read(10,*,end=22)az1,ip1
      if(az1(5).le.0.001)zz=az1(3)
      if(az1(5).gt.0.001)zz=az1(5)
      
      zmag0=0.0
      do i=1,N0-1
         if(zz.ge.z3(i).and.zz.lt.z3(i+1))then
            zmag0=z4(i)+(z4(i+1)-z4(i))*(zz-z3(i))/(z3(i+1)-z3(i))
         endif
      enddo
      
      zmass=10**(az1(18)-10.)
      if(zmass.ge.0.8)then
         N=N+1
         z1(N,1)=az1(1)
         z1(N,2)=az1(2)
         if(az1(5).le.0.001)then
            index1(N)=0
            z1(N,3)=az1(3)
         endif
         if(az1(5).gt.0.001)then
            index1(N)=1
            z1(N,3)=az1(5)
         endif
         z1(N,4)=az1(12)
         z1(N,5)=az1(14)
         z1(N,6)=zmass
         z1(N,7)=zmag0
      endif
      goto 11
 22   close(10)
      write(*,*)N

 33   read(20,*,err=33,end=44)id2,az2,ip2
      zmstar=0.0
      do i=1,N0-1
         if(az2(3).ge.z3(i).and.az2(3).lt.z3(i+1))then
            zmstar=z4(i)+(z4(i+1)-z4(i))*(az2(3)-z3(i))/(z3(i+1)-z3(i))
         endif
      enddo
      if(az2(3).le.0.7)zgap=0.04*(1+az2(3))
      if(az2(3).gt.0.7)zgap=0.15*az2(3)-0.037
      ddz=dis(az2(3))
      
      sum=0.
      M=0
      M2=0
      do i=1,N
         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.1
     $        .and.ip2(2).eq.1.and.
     $        abs(az2(3)-z1(i,3)).le.0.00833*(1+az2(3)))then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=ddz*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zmstar-z1(i,7)))
            if(z6.le.az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M=M+1
               sum=sum+zmass2
            endif
            if(z6.le.0.5*az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M2=M2+1
            endif
         endif

         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.1
     $        .and.ip2(2).eq.0.and.
     $        abs(az2(3)-z1(i,3)).le.zgap)then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=ddz*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zmstar-z1(i,7)))
            if(z6.le.az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M=M+1
               sum=sum+zmass2
            endif
            if(z6.le.0.5*az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M2=M2+1
            endif
         endif

         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.0
     $        .and.abs(az2(3)-z1(i,3)).le.zgap)then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=ddz*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zmstar-z1(i,7)))
            if(z6.le.az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M=M+1
               sum=sum+zmass2
            endif
            if(z6.le.0.5*az2(7).and.zmass2.ge.1..and.
     $           z1(i,4).ge.az2(4)-0.25)then
               M2=M2+1
            endif
         endif
         
      enddo
      rich=sum-az2(10)*az2(7)**2
      if(rich.ge.20..and.M.ge.5)then
         write(30,100)id2,(az2(j),j=1,9),rich,ip2(1),M,ip2(2)
      endif
      goto 33
 44   close(20)

      goto 15
 16   close(12)
 100  format(I7,2F11.5,F8.4,2F7.3,F7.2,F6.3,3F8.2,3I4)
      end

c     id,ra,dec,z,i_bcg,w1_bcg,stamass_bcg,r500,lum05e,lum05,lum500,N5e,N500
      
      function dis(z)
      real dis,z
      N_step=100
      dz=z/N_step
      sum=0.0
      do j=1,N_step
         sum=sum+dz/sqrt(0.7+0.3*(1+dz*j)**3)
      enddo
      dis=4285.7*sum/(1+z)/57.3 !h=0.7
      end
