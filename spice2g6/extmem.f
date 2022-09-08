      subroutine extmem(ipntr,ksize)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
c
c***  extmem - extend size of existing block
c
c
c...  check for valid pointer
      if (memptr(ipntr(1))) go to 10
      memerr=5
      call errmem(2,memerr,ipntr(1))
   10 isize=ksize*istack(ltab+5)
c...  check for valid size
      if (isize.ge.0) go to 20
      memerr=2
      call errmem(2,memerr,ipntr(1))
c...  check if enough space already there
   20 if ((istack(ltab+2)-istack(ltab+3)).ge.isize) go to 40
      need=nxtevn(isize)-memavl
      if (need.le.0) go to 30
c...  insufficient space -- bump memory size
      need=nxtmem(need)
      icore=icore+need
      call memory
      if(memerr.ne.0) call errmem(2,memerr,ipntr(1))
      ltab1=ldval-ntab
      istack(ltab1+2)=istack(ltab1+2)+need
c...  relocate block entry table
      nwords=numblk*ntab
      cpyknt=cpyknt+dble(nwords)
      call copy4(istack(loctab+1),istack(loctab+need+1),nwords)
      loctab=loctab+need
      ldval=ldval+need
      memavl=memavl+need
      ltab=ltab+need
c...  move blocks to make space
   30 continue
      call comprs(0,ltab)
      call comprs(1,ltab)
   40 jsize=istack(ltab+3)
      istack(ltab+3)=istack(ltab+3)+isize
      memavl=memavl-(nxtevn(istack(ltab+3))-nxtevn(jsize))
      call memadj
      return
      end
