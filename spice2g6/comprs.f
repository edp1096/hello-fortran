      subroutine comprs(icode,limit)
      implicit double precision (a-h,o-z)
c
c      this routine compresses all available memory into a single block.
c if *icode* is zero, compression of memory from word 1 to *limit* is
c done;  otherwise, compression from *ldval* down to *limit* is done.
c
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      if (icode.ne.0) go to 100
      nblk=numblk
      ltab2=loctab
   10 ltab1=ltab2
      if (ltab1.ge.limit) go to 200
      if (nblk.eq.1) go to 200
      nblk=nblk-1
      ltab2=ltab1+ntab
      morg=istack(ltab1+1)
      msiz=istack(ltab1+2)
      muse=nxtevn(istack(ltab1+3))
      mslp=istack(ltab1+6)
      if ((msiz-muse).le.mslp) go to 10
      muse=muse+mslp
c...  move succeeding block down
      morg2=istack(ltab2+1)
      muse2=istack(ltab2+3)
      madr2=istack(ltab2+4)
      iwsize=istack(ltab2+5)
      if (madr2.ne.0) go to 15
      if (muse2.eq.0) go to 20
   15 cpyknt=cpyknt+dble(muse2)
      call copy4(istack(nwoff+morg2+1),istack(nwoff+morg+muse+1),muse2)
      istack(lorg+madr2)=(morg+muse)/iwsize
   20 istack(ltab1+2)=muse
      istack(ltab2+1)=morg+muse
      istack(ltab2+2)=istack(ltab2+2)+(msiz-muse)
      go to 10
c
c
  100 nblk=numblk
      ltab2=ldval-ntab
  110 ltab1=ltab2
      if (ltab1.le.limit) go to 200
      if (nblk.eq.1) go to 200
      nblk=nblk-1
      ltab2=ltab1-ntab
      morg=istack(ltab1+1)
      msiz=istack(ltab1+2)
      muser=istack(ltab1+3)
      muse=nxtevn(muser)
      madr=istack(ltab1+4)
      iwsize=istack(ltab1+5)
      mslp=istack(ltab1+6)
      if ((msiz-muse).le.mslp) go to 110
      muse=muse+mslp
      mspc=msiz-muse
      cpyknt=cpyknt+dble(muser)
      call copy4(istack(nwoff+morg+1),istack(nwoff+morg+mspc+1),muser)
      istack(ltab1+1)=morg+mspc
      istack(ltab1+2)=muse
      istack(ltab2+2)=istack(ltab2+2)+mspc
      if (madr.eq.0) go to 110
      istack(lorg+madr)=(morg+mspc)/iwsize
      go to 110
c...  all done
  200 return
      end
