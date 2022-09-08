      logical function memptr(ipntr)
      implicit double precision (a-h,o-z)
c
c      this routine checks whether *ipntr* is a valid block pointer.
c if it is valid, *ltab* is set to point to the corresponding entry in
c the block table.
c
c... ipntr is an array to avoid 'call by value' problems (see setmem)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      memptr=.false.
      ltab=loctab
      locpnt=locf(ipntr(1))
      do 20 i=1,numblk
      if (locpnt.ne.istack(ltab+4)) go to 10
      if (ipntr(1)*istack(ltab+5).ne.istack(ltab+1)) go to 10
      memptr=.true.
      go to 30
   10 ltab=ltab+ntab
   20 continue
   30 return
      end
