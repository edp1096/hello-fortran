      subroutine dmpmem(ipntr)
      implicit double precision (a-h,o-z)
c
c      this routine prints out the current memory allocation map.
c *ipntr* is the table pointer of the current memory manager call
c
c spice version 2g.6  sccsid=tabinf 3/15/83
      common /tabinf/ ielmnt,isbckt,nsbckt,iunsat,nunsat,itemps,numtem,
     1   isens,nsens,ifour,nfour,ifield,icode,idelim,icolum,insize,
     2   junode,lsbkpt,numbkp,iorder,jmnode,iur,iuc,ilc,ilr,numoff,isr,
     3   nmoffc,iseq,iseq1,neqn,nodevs,ndiag,iswap,iequa,macins,lvnim1,
     4   lx0,lvn,lynl,lyu,lyl,lx1,lx2,lx3,lx4,lx5,lx6,lx7,ld0,ld1,ltd,
     5   imynl,imvn,lcvn,nsnod,nsmat,nsval,icnod,icmat,icval,
     6   loutpt,lpol,lzer,irswpf,irswpr,icswpf,icswpr,irpt,jcpt,
     7   irowno,jcolno,nttbr,nttar,lvntmp
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      dimension ipntr(1)
      dimension aptr(75)
      data aptr /6hielmnt,6hisbckt,6hnsbckt,6hiunsat,6hnunsat,6hitemps,
     1 6hnumtem,6hisens ,6hnsens ,6hifour ,6hnfour ,6hifield,
     2 6hicode ,6hidelim,6hicolum,6hinsize,
     3 6hjunode,6hlsbkpt,6hnumbkp,6hiorder,6hjmnode,
     4 6hiur   ,6hiuc   ,6hilc   ,6hilr   ,6hnumoff,6hisr   ,
     5 6hnmoffc,6hiseq  ,6hiseq1  ,6hneqn  ,6hnodevs,
     6 6hndiag ,6hiswap ,6hiequa ,6hmacins,6hlvnim1,
     7 6hlx0   ,6hlvn   ,6hlynl  ,6hlyu   ,6hlyl   ,
     8 6hlx1   ,6hlx2   ,6hlx3   ,6hlx4   ,6hlx5   ,6hlx6   ,
     9 6hlx7   ,6hld0   ,6hld1   ,6hltd   ,6himynl ,6himvn  ,6hlcvn  ,
     * 6hnsnod ,6hnsmat ,6hnsval ,6hicnod ,6hicmat ,6hicval ,6hloutpt,
     * 6hlpol  ,6hlzer  ,6hirswpf,6hirswpr,6hicswpf,6hicswpr,6hirpt  ,
     * 6hjcpt  ,6hirowno,6hjcolno,6hnttbr ,6hnttar ,6hlvntmp/
      data ablnk /1h /
c
      iaddr=locf(ielmnt)-1
      itemp=locf(ipntr(1))-iaddr
      anam=ablnk
      if(itemp.gt.0.and.itemp.le.75) anam=aptr(itemp)
      iadr=locf(ipntr(1))
      write (iofile,5) anam,iadr,icore,maxmem,memavl,ldval
    5 format('0current pointer ',a6,'@ = z',i6,/' corsiz=',i7,
     1  /' maxmem=',i7,/' avlspc=',i7,/' ldval=',i7,
     2  /1h0,24x,'memory allocation map'/14x,'blknum memorg memsiz',
     3  '  memuse usrptr  addr    name')
      ltab1=loctab
      do 20 i=1,numblk
      morg=istack(ltab1+1)
      msiz=istack(ltab1+2)
      muse=istack(ltab1+3)
      madr=istack(ltab1+4)
      anam=ablnk
      ndex=madr-iaddr
      if(ndex.gt.0.and.ndex.le.75) anam=aptr(ndex)
      jptr=0
      if (madr.gt.0) jptr=istack(lorg+madr)
      write (iofile,11) i,morg,msiz,muse,jptr,madr,anam
   11 format(13x,5i7,3x,i7,'z',1x,a6)
      ltab1=ltab1+ntab
   20 continue
      write (iofile,21)
   21 format(1h0,24x,'end of allocation map'/)
      return
      end
