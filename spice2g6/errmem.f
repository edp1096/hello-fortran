      subroutine errmem(inam,ierror,ipntr)
      implicit double precision (a-h,o-z)
      dimension ipntr(1)
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      dimension errnam(7)
      data errnam /6hclrmem,6hextmem,6hgetmem,6hptrmem,6hrelmem,
     1   6hsetmem,6hsizmem/
c
      go to (200,410,420,300,510,530),ierror
c
c*** error(s) found ***
c
c.. nxtevn and/or nxtmem incompatible with nwd4, nwd8, and nwd16
c
  200 write(iofile,201)
  201 format('0memory manager variables nwd4-8-16 incompatible with nxte
     1vn and nxtmem')
      go to 900
c
c...  memory needs exceed maximum available space
  300 write (iofile,301) maxmem
  301 format('0*error*:  memory requirement exceeds machine capacity',
     1/'0 memory needs exceed',i6)
      go to 900
c...    *isize* < 0
  410 write(iofile,411)
  411 format('0size parameter negative')
      go to 900
c...  getmem:  attempt to reallocate existing block
  420 write(iofile,421)
  421 format('0attempt to reallocate existing table')
      go to 900
c...    *ipntr* invalid
  510 write(iofile,511)
  511 format('0table pointer invalid')
      go to 900
c...  relmem:  *isize* larger than indicated block
  530 write(iofile,531)
  531 format('0attempt to release more than total table')
c...  issue error message
  900 write (iofile,901) errnam(inam)
  901 format('0*abort*:  internal memory manager error at entry ',
     1  a7)
  950 call dmpmem(ipntr(1))
 1000 stop
      end
