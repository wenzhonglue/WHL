      parameter (ang=3.14159265/180.)
      character fn1*18,fn2*18,aa*100
      real az1(18),z1(8000000,5),z2(8000000,5)
      real a5(12),z3(200),z4(200)
      integer id1,ip1,index2(8000000)
      open(12,file='cand_bcg.lst')
 15   read(12,*,end=16)fn1,fn2
      open(10,file='../data/'//fn1)
      open(20,file=fn2)
      N=0
      M=0
 11   read(10,*,end=22)az1,ip1
      
      if(az1(18).ge.10.255)then
         N=N+1
         z1(N,1)=az1(1)
         z1(N,2)=az1(2)
         z1(N,3)=az1(3)
         z1(N,4)=az1(12)
         z1(N,5)=az1(18)
      endif
      if(ip1.eq.1)then
         M=M+1
         z2(M,1)=az1(1)
         z2(M,2)=az1(2)
         if(az1(5).le.0.001)then
            index2(M)=0
            z2(M,3)=az1(3)
         endif
         if(az1(5).gt.0.001)then
            index2(M)=1
            z2(M,3)=az1(5)
         endif         
         z2(M,4)=az1(12)
         z2(M,5)=az1(18)
      endif
      goto 11
 22   close(10)
      write(*,*)N,M

      do i=1,M
         n50=0
         if(z2(i,3).le.0.7)zgap=0.04*(1+z2(i,3))
         if(z2(i,3).gt.0.7)zgap=0.15*z2(i,3)-0.037
         do j=1,N
            if(abs(z2(i,1)-z1(j,1)).le.0.5.and.abs(z2(i,2)-z1(j,2))
     $           .le.0.5.and.abs(z2(i,3)-z1(j,3)).le.zgap !)then
     $           .and.z1(j,5).le.z2(i,5)+0.1)then
               z7=abs(z2(i,1)-z1(j,1))**2*cos(z2(i,2)*ang)**2+
     $              abs(z2(i,2)-z1(j,2))**2
               z6=dis(z2(i,3))*sqrt(z7)
               if(z6.le.0.5)n50=n50+1
            endif
         enddo
         if(n50.ge.3)then
            write(20,100)i,(z2(i,k),k=1,5),n50,index2(i)
         endif
      enddo

      goto 15
 16   close(12)
 100  format(I6,2F11.5,F8.4,2F7.3,2I4)
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
