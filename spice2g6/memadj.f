      subroutine memadj
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
   50 maxuse=max0(maxuse,(ldval-memavl-ifwa))
      memdec=2*nxtmem(1)
      if (memavl.lt.memdec) return
c...  compress current allocations of memory
      call comprs(0,ldval)
c...  adjust memory size
      memdel=0
   60 icore=icore-memdec
      memdel=memdel+memdec
      memavl=memavl-memdec
      if (memavl.ge.memdec) go to 60
      ltab1=ldval-ntab
      istack(ltab1+2)=istack(ltab1+2)-memdel
c...  relocate block entry table
      nwords=numblk*ntab
      cpyknt=cpyknt+dble(nwords)
      call copy4(istack(loctab+1),istack(loctab-memdel+1),nwords)
      loctab=loctab-memdel
      ldval=ldval-memdel
      call memory
      return
      end
