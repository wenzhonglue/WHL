      parameter (ang=3.14159265/180.)
      character aa*100
      character fn1*18,fn2*18,fn3*18
      real az1(18),az2(5),z1(8000000,7),tmp(200,6)
      real az4(2),z2(50000,2)
      real a5(12),z3(200),z4(200),a5x(12),z3x(200),z4x(200)
      real a6(10),y3(200),y4(200)
      integer id1,id2,ip1,ip2(2),id4,index1(8000000)
      open(12,file='clusterbcg.lst')
      open(25,file='bc2003_magnitude_F119')
      open(40,file='../data/visual.dat')
      N1=0
      do i=1,50
         read(25,'(a)')aa
      enddo
 26   read(25,*,end=27)a5x
      N1=N1+1
      z3x(N1)=a5x(1)
      z4x(N1)=a5x(12)-28.
      goto 26
 27   close(25)
      write(*,*)N1

      N3=0
 55   read(40,*,err=55,end=66)id4,az4
      N3=N3+1
      z2(N3,1)=az4(1)
      z2(N3,2)=az4(2)
      goto 55
 66   close(40)
      write(*,*)N3
      
 15   read(12,*,end=16)fn1,fn2,fn3
      open(10,file='../data/'//fn1)
      open(20,file=fn2)
      open(30,file=fn3)
      N=0
 11   read(10,*,end=22)az1,ip1
      if(az1(5).le.0.001)zz=az1(3)
      if(az1(5).gt.0.001)zz=az1(5)
      
      zmag0=0.
      do i=1,N1-1
         if(zz.ge.z3x(i).and.zz.lt.z3x(i+1))then
            zmag0=z4x(i)+
     $           (z4x(i+1)-z4x(i))*(zz-z3x(i))/(z3x(i+1)-z3x(i))
         endif
      enddo

c     --------------------------------------------
      do i=1,N3
         if(abs(az1(1)-z2(i,1)).le.0.01.and.abs(az1(2)-z2(i,2))
     $        .le.0.01)then
            z9=abs(az1(1)-z2(i,1))**2*cos(az1(2)*ang)**2
     $           +abs(az1(2)-z2(i,2))**2
            z8=3600.0*sqrt(z9)
            if(z8.le.2.)ip1=0
         endif
      enddo
      
      if(ip1.eq.1)then
         zmass=10**(az1(18)-10.)
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

 33   read(20,*,end=44)id2,az2,ip2
c     ----------------------------
c      if(az2(3).ge.0.9.and.az2(3).le.1.15)az2(3)=0.8*az2(3)+0.18
c     ----------------------------
      zlim=0.0
      do i=1,N1-1
         if(az2(3).ge.z3x(i).and.az2(3).lt.z3x(i+1))then
            zlim=z4x(i)+
     $           (z4x(i+1)-z4x(i))*(az2(3)-z3x(i))/(z3x(i+1)-z3x(i))
         endif
      enddo
      if(az2(3).le.0.7)zgap=0.04*(1+az2(3))
      if(az2(3).gt.0.7)zgap=0.15*az2(3)-0.037
      rad=0.5
      
      M=0
      do i=1,N
         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.1
     $        .and.ip2(2).eq.1.and.
     $        abs(az2(3)-z1(i,3)).le.0.00833*(1+az2(3)))then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=dis(az2(3))*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zlim-z1(i,7)))
            if(z6.le.rad.and.zmass2.gt.10.)then
               M=M+1
               tmp(M,1)=z1(i,1)
               tmp(M,2)=z1(i,2)
               tmp(M,3)=z1(i,3)
               tmp(M,4)=z1(i,4)
               tmp(M,5)=z1(i,5)
               tmp(M,6)=zmass2
            endif
         endif

         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.1
     $        .and.ip2(2).eq.0.and.
     $        abs(az2(3)-z1(i,3)).le.0.025*(1+az2(3)))then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=dis(az2(3))*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zlim-z1(i,7)))
            if(z6.le.rad.and.zmass2.gt.10.)then
               M=M+1
               tmp(M,1)=z1(i,1)
               tmp(M,2)=z1(i,2)
               tmp(M,3)=z1(i,3)
               tmp(M,4)=z1(i,4)
               tmp(M,5)=z1(i,5)
               tmp(M,6)=zmass2
            endif
         endif

         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.0.and.
     $        abs(az2(3)-z1(i,3)).le.0.6*zgap)then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=dis(az2(3))*sqrt(z7)
            zmass2=z1(i,6)*10**(0.4*(zlim-z1(i,7)))
            if(z6.le.rad.and.zmass2.gt.10.)then
               M=M+1
               tmp(M,1)=z1(i,1)
               tmp(M,2)=z1(i,2)
               tmp(M,3)=z1(i,3)
               tmp(M,4)=z1(i,4)
               tmp(M,5)=z1(i,5)
               tmp(M,6)=zmass2
            endif
         endif
      enddo
      if(M.ge.1)then
         chi2=99.
         do j=1,M
            if(tmp(j,4).le.chi2)then
               chi2=tmp(j,4)
               k=j
            endif
         enddo
         write(30,100)id2,tmp(k,1),tmp(k,2),az2(3),tmp(k,4),tmp(k,5),
     $        tmp(k,6),ip2
      endif
      goto 33
 44   close(20)
      
      goto 15
 16   close(12)
 100  format(I7,2F11.5,F8.4,2F7.3,F7.2,2I4)
      end

c     id,ra,dec,z,z_bcg,w1_bcg,stamass_bcg
      
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
      
      function ez(z)
      real ez,z
      ez=0.7+0.3*(1+z)**3
      end
      
