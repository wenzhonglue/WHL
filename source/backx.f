      parameter (ang=3.14159265/180.)
      character fn1*18,fn2*18,fn3*18,aa*100
      real az1(18),az2(6),z1(8000000,7),tmp(500,2),z2(500),err1(500)
      real a5(12),z3(200),z4(200)
      integer id1,id2,ip1,ip2(2),index1(8000000),ng2(500000)
      open(12,file='back.lst')
 15   read(12,*,end=16)fn1,fn2,fn3
      open(10,file='../data/'//fn1)
      open(20,file=fn2)
      open(30,file=fn3)
      N=0
 11   read(10,*,end=22)az1,ip1
      
      if(az1(18).ge.10.)then
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
         z1(N,6)=10**(az1(18)-10.)
         z1(N,7)=az1(4)
         if(z1(N,7).le.0.005)z1(N,7)=0.005
      endif
      goto 11
 22   close(10)
      write(*,*)N

 33   read(20,*,err=33,end=44)id2,az2,ip2
      if(az2(3).le.0.7)zgap=0.04*(1+az2(3))
      if(az2(3).gt.0.7)zgap=0.15*az2(3)-0.037
      rad=1.0/ez(az2(3))**0.333
      ddz=dis(az2(3))
      M=0
      do i=1,N
         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $        abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.0
     $        .and.abs(az2(3)-z1(i,3)).le.zgap)then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=ddz*sqrt(z7)
            if(z6.le.rad.and.z1(i,4).ge.az2(4)-0.25)then
               M=M+1
               tmp(M,1)=z1(i,4)
               tmp(M,2)=z1(i,6)
            endif
         endif
      enddo
      zw2=99.
      zm2=0.
      if(M.ge.1)then
         do j=1,M
            if(tmp(j,1).ge.az2(4)+0.05.and.tmp(j,1).lt.zw2)then
               zw2=tmp(j,1)
               zm2=tmp(j,2)
            endif
         enddo
      endif
      
      sum=0.
      do i=1,N
         if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.1..and.
     $        abs(az2(2)-z1(i,2)).le.1..and.
     $        abs(az2(3)-z1(i,3)).le.zgap)then
            z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $           abs(az2(2)-z1(i,2))**2
            z6=ddz*sqrt(z7)
            if(z6.ge.2..and.z6.le.4..and.z1(i,4).ge.zw2)then
               sum=sum+z1(i,6)
            endif
         endif
      enddo
      flat=sum/12.
      write(30,100)id2,az2,zw2,zm2,flat,ip2
      goto 33
 44   close(20)

      goto 15
 16   close(12)
 100  format(I6,2F11.5,F8.4,2F7.3,F7.2,F7.3,2F7.2,2I4)
      end
      
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
      
