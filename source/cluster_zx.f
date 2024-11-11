      parameter (ang=3.14159265/180.)
      character fn1*18,fn2*18,fn3*18,aa*100
      real az1(18),az2(5),z1(8000000,5),z2(500000,5),tmp(500),tmp2(500)
      real a5(12),z3(200),z4(200)
      integer id1,id2,ip1,ip2(2),index1(8000000),ng2(500000)
      open(12,file='cand_bcg.lst')
 15   read(12,*,end=16)fn1,fn2,fn3
      open(10,file='../data/'//fn1)
      open(20,file=fn2)
      open(30,file=fn3)
      N=0
 11   read(10,*,end=22)az1,ip1
      
      if(az1(18).ge.10.75)then
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
         
         z1(N,4)=az1(4)
         z1(N,5)=az1(12)

      endif
      goto 11
 22   close(10)
      write(*,*)N

      M=0
 33   read(20,*,end=44)id2,az2,ip2
      if(az2(3).le.0.7)zgap=0.04*(1+az2(3))
      if(az2(3).gt.0.7)zgap=0.15*az2(3)-0.037
      rad=1.0/ez(az2(3))**0.333
      ddz=dis(az2(3))
      zz=0.

      if(ip2(2).eq.1)zz=az2(3)

      if(ip2(2).eq.0)then
         k1=0
         do i=1,N
            if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $           abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.1
     $           .and.abs(az2(3)-z1(i,3)).le.0.025*(1+az2(3)))then
               z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $              abs(az2(2)-z1(i,2))**2
               z6=ddz*sqrt(z7)
               if(z6.le.rad)then
                  k1=k1+1
                  tmp(k1)=z1(i,3)
               endif
            endif
         enddo
         k2=0
         do i=1,N
            if(abs(az2(1)-z1(i,1))*cos(az2(2)*ang).le.0.5.and.
     $           abs(az2(2)-z1(i,2)).le.0.5.and.index1(i).eq.0
     $           .and.abs(az2(3)-z1(i,3)).le.0.6*zgap)then
               z7=abs(az2(1)-z1(i,1))**2*cos(az2(2)*ang)**2+
     $              abs(az2(2)-z1(i,2))**2
               z6=ddz*sqrt(z7)
               if(z6.le.rad)then
                  k2=k2+1
                  tmp2(k2)=z1(i,3)
               endif
            endif
         enddo

         if(k1.ge.1)then
            ip2(2)=1
            M=M+1
            if(mod(k1,2).eq.1)zz=select(k1/2+1,k1,tmp)
            if(mod(k1,2).eq.0)then
               zz=0.5*(select(k1/2,k1,tmp)+select(k1/2+1,k1,tmp))
            endif
         endif
         if(k1.eq.0.and.k2.ge.1)then
            if(mod(k2,2).eq.1)zz=select(k2/2+1,k2,tmp2)
            if(mod(k2,2).eq.0)then
               zz=0.5*(select(k2/2,k2,tmp2)+select(k2/2+1,k2,tmp2))
            endif
         endif
      endif

      if(zz.ge.0.001)then
         write(30,100)id2,az2(1),az2(2),zz,az2(4),az2(5),ip2
      endif
      goto 33
 44   close(20)
      write(*,*)M
      
      goto 15
 16   close(12)
 100  format(I7,2F11.5,F8.4,2F7.3,2I4)
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
      
      include "select.for"
