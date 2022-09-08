      subroutine relmem(ipntr,ksize)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
c
c***  relmem - release part of block
c
c
c...  check for valid pointer
      if (memptr(ipntr(1))) go to 10
      memerr=5
      call errmem(5,memerr,ipntr(1))
   10 isize=ksize*istack(ltab+5)
c...  check for valid size
      if (isize.ge.0) go to 20
      memerr=2
      call errmem(5,memerr,ipntr(1))
   20 jsize=istack(ltab+3)
      if (isize.le.jsize) go to 30
      memerr=6
      call errmem(5,memerr,ipntr(1))
   30 istack(ltab+3)=istack(ltab+3)-isize
      memavl=memavl+(nxtevn(jsize)-nxtevn(istack(ltab+3)))
      call memadj
      return
      end
