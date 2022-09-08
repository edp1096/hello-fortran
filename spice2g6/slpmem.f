      subroutine slpmem(ipntr,ksize)
      implicit double precision (a-h,o-z)
c
c      this routine may be used to define a certain amount of `slop' to
c be associated with a particular table managed by the memory manager.
c this *slop* is defined as a number of entries in the table for which
c space is to be held ***if possible*** during compaction of the managed
c area of memory.  this feature can eliminate the overhead incurred by
c alternatively extending more than one table at a time.  (for example,
c if the program contains a code sequence
c
c                  do 100 i=1,500
c                     ...
c                  call extmem(table1,1)
c                     ...
c                  call extmem(table2,1)
c                     ...
c              100 continue
c
c then the overhead incurred by this memory manager can be reduced to
c essentially nothing if prior to the above code sequence the program
c executes
c
c                  call slpmem(table1,20)
c                  call slpmem(table2,20)
c
c where `20' is a typical number (for the above example, the memory-to-
c memory copying overhead of the memory manager would be reduced by a
c factor of 20).
c
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
      dimension ipntr(1)
      islp=nxtevn(ksize)
      call extmem(ipntr,islp)
      islp=islp*istack(ltab+5)
      istack(ltab+3)=istack(ltab+3)-islp
      memavl=memavl+istack(ltab+6)
      istack(ltab+6)=islp
      return
      end
