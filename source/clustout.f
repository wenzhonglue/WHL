      character fn1*18,fn2*18,fn3*18
      real az1(5),az2(10),z1(100,5)
      integer id2,ip2(3)
      open(30,file='clustout.tmp')
      open(12,file='cluslumx.lst')
 15   read(12,*,end=16)fn1,fn2,fn3
      open(20,file=fn3)
 33   read(20,*,end=44)id2,az2,ip2

      if(ip2(2).ge.5)then
         write(30,100)id2,(az2(j),j=1,7),az2(10),ip2(2),ip2(3)
      endif
      goto 33
 44   close(20)

      goto 15
 16   close(12)
 100  format(I7,2F11.5,F8.4,2F7.3,F7.2,F6.3,F8.2,2I4)
      end
      
c     id,ra,dec,z,i_bcg,w1_bcg,stamass_bcg,r500,lum500,Nr500
