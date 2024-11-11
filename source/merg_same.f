      parameter (ang=3.14159265/180.)
      real az1(8),z1(8000000,8)
      integer id1,ip1(2),ip2(8000000,2),index1(8000000)
      open(10,file='clustout.tmp')
      open(20,file='cluster.dat')
      N=0
 11   read(10,*,end=22)id1,az1,ip1
      rich=az1(8)*(1+az1(3))**0.201
      if(rich.ge.40..and.ip1(1).ge.6)then
         N=N+1
         do j=1,7
            z1(N,j)=az1(j)
         enddo
         z1(N,8)=rich
         ip2(N,1)=ip1(1)
         ip2(N,2)=ip1(2)
         index1(N)=1
      endif
      goto 11
 22   close(10)
      write(*,*)N

      do i=1,N
         do j=1,N
            zx=0.5*(z1(i,3)+z1(j,3))
            if(zx.le.0.7)zgap=0.04*(1+zx)
            if(zx.gt.0.7)zgap=0.15*zx-0.037
            if(z1(i,7).ge.z1(j,7))rij=z1(i,7)
            if(z1(i,7).le.z1(j,7))rij=z1(j,7)
            if(i.lt.j.and.abs(z1(i,3)-z1(j,3)).le.1.5*zgap
     $           .and.abs(z1(i,1)-z1(j,1)).le.0.5
     $           .and.abs(z1(i,2)-z1(j,2)).le.0.5)then
               z7=abs(z1(i,1)-z1(j,1))**2*cos(z1(i,2)*ang)**2+
     $              abs(z1(i,2)-z1(j,2))**2
               z6=dis(zx)*sqrt(z7)
               if(z6.le.1.5*rij)then
                  if(z1(i,8).gt.z1(j,8))index1(j)=0
                  if(z1(i,8).le.z1(j,8))index1(i)=0
               endif
            endif
         enddo
      enddo
      
      M=0
      do i=1,N
         if(index1(i).eq.1)then
            M=M+1
            write(20,100)M,(z1(i,j),j=1,8),ip2(i,1),ip2(i,2)
         endif
      enddo
 100  format(I7,2F11.5,F8.4,2F7.3,F7.2,F6.3,F8.2,2I4)
      end

c     id,ra,dec,z,i_bcg,w1_bcg,stamass_bcg,r500,lum500,Nr500
      
      function dis(z)
      real dis,z
      N_step=100
      dz=z/N_step
      sum=0.0
      do j=1,N_step
         sum=sum+dz/sqrt(0.7+0.3*(1+dz*j)**3)
      enddo
      dis=4285.7*sum/(1+z)/57.3
      end

      function dmode(z)
      real z,ldis
      n=1000
      dz=z/real(n)
      sum=0.
      do i=1,n
         sum=sum+dz/sqrt(0.7+0.3*(1+i*dz)**3)
      enddo
      tmp=4285.7*sum*(1+z)
      dmode=25.+5.*log10(tmp)
      end
