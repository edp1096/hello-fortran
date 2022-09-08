      subroutine ptrmem(ipntr,ipntr2)
      implicit double precision (a-h,o-z)
      dimension ipntr(1),ipntr2(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
c
c***  ptrmem - reset memory pointer
c
c...  verify that pointer is valid
      if (memptr(ipntr(1))) go to 10
      memerr=5
      call errmem(4,memerr,ipntr(1))
c...  reset block pointer to be *ipntr2*
   10 ipntr2(1)=ipntr(1)
      istack(ltab+4)=locf(ipntr2(1))
      call memadj
      return
      end
