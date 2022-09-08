      subroutine cpytb8(itabo,itabn)
      implicit double precision (a-h,o-z)
c
c     this routine copies a table.  its use is made necessary by the
c fact that only one pointer is allowed per table.
c
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      call sizmem(nodplc(itabo),isize)
      call getm8(nodplc(itabn),isize)
      loco=nodplc(itabo)
      locn=nodplc(itabn)
      call copy8(value(loco+1),value(locn+1),isize)
      return
      end
