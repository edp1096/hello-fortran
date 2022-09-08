      subroutine sizmem(ipntr,ksize)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
c
c***  sizmem - determine size of existing block
c
c
c...  check for valid pointer
      if (memptr(ipntr(1))) go to 10
      memerr=5
      call errmem(7,memerr,ipntr(1))
   10 ksize=istack(ltab+3)/istack(ltab+5)
      return
      end
