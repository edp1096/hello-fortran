      subroutine clrmem(ipntr)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
c
c***  clrmem - release block
c
c
c...  check that pointer is valid
      if (memptr(ipntr(1))) go to 10
      memerr=5
      call errmem(1,memerr,ipntr(1))
   10 msiz=istack(ltab+2)
      muse=istack(ltab+3)
      memavl=memavl+nxtevn(muse)+istack(ltab+6)
c...  assumption:  first allocated block is never cleared.
      ltab1=ltab-ntab
      istack(ltab1+2)=istack(ltab1+2)+msiz
c...  reposition the block table
      nwords=ltab-loctab
      cpyknt=cpyknt+dble(nwords)
      call copy4(istack(loctab+1),istack(loctab+ntab+1),nwords)
      numblk=numblk-1
      loctab=loctab+ntab
      memavl=memavl+ntab
      ltab1=ldval-ntab
      istack(ltab1+2)=istack(ltab1+2)+ntab
      ipntr(1)=2**30-1
      call memadj
      return
      end
