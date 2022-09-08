      subroutine getmx(ipntr,ksize,iwsize)
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      logical memptr
      dimension ipntr(1)
c
c***  getmem - get block
c
c
      isize=ksize*iwsize
c...  check for valid size
      if (isize.ge.0) go to 5
      memerr=2
      call errmem(3,memerr,ipntr(1))
c...  check for attempt to reallocate existing block
    5 if (.not.memptr(ipntr(1))) go to 8
      memerr=3
      call errmem(3,memerr,ipntr(1))
    8 jsize=nxtevn(isize)
      call comprs(0,ldval)
c...  check if enough space already there
      need=jsize+ntab-memavl
      if (need.le.0) go to 10
c...  insufficient space -- bump memory size
      need=nxtmem(need)
      icore=icore+need
      call memory
      if(memerr.ne.0) call errmem(3,memerr,ipntr(1))
      ltab1=ldval-ntab
      istack(ltab1+2)=istack(ltab1+2)+need
c...  relocate block entry table
      nwords=numblk*ntab
      cpyknt=cpyknt+dble(nwords)
      call copy4(istack(loctab+1),istack(loctab+need+1),nwords)
      loctab=loctab+need
      ldval=ldval+need
      memavl=memavl+need
c...  a block large enough now exists -- allocate it
   10 ltab1=ldval-ntab
      morg=istack(ltab1+1)
      msiz=istack(ltab1+2)
      muse=istack(ltab1+3)
      muse=nxtevn(muse)
      madr=istack(ltab1+4)
c...  construct new table entry
   15 istack(ltab1+2)=muse
      loctab=loctab-ntab
      nwords=numblk*ntab
      cpyknt=cpyknt+dble(nwords)
      call copy4(istack(loctab+ntab+1),istack(loctab+1),nwords)
      numblk=numblk+1
      memavl=memavl-ntab
      istack(ltab1+1)=morg+muse
      istack(ltab1+2)=msiz-muse-ntab
c...  set user size into table entry for this block
   20 istack(ltab1+3)=isize
      istack(ltab1+4)=locf(ipntr(1))
      istack(ltab1+5)=iwsize
      istack(ltab1+6)=0
      memavl=memavl-jsize
      ipntr(1)=istack(ltab1+1)/iwsize
      call memadj
      return
      end
