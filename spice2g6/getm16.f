      subroutine getm16(ipntr,ksize)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      iwsize=nwd16
      call getmx(ipntr(1),ksize,iwsize)
      return
      end
