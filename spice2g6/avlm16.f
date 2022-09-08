      subroutine avlm16(iavl)
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      iavl=((maxmem-icore)/nxtmem(1))*nxtmem(1)-ntab+memavl
      iavl=iavl/nwd16
      return
      end
