      subroutine setplt(loc)
      implicit double precision (a-h,o-z)
c
c     this routine generates the 'legend' subheading used to identify
c individual traces on multi-trace line-printer plots.
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=outinf 3/15/83
      common /outinf/ xincr,string(15),xstart,yvar(8),itab(8),itype(8),
     1   ilogy(8),npoint,numout,kntr,numdgt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      dimension logopt(6)
      data logopt / 2, 2, 1, 1, 1, 1 /
      data ablnk, atimex, afreq / 1h , 6h  time, 6h  freq /
      data pltsym / 8h*+=$0<>? /
c
c  set limits depending upon the analysis mode
c
      if (mode-2) 10,20,30
   10 xstart=tcstar(1)
      xincr=tcincr(1)
      npoint=icvflg
      itemp=itcelm(1)
      loce=nodplc(itemp+1)
      asweep=value(loce)
      go to 40
   20 xstart=tstart
      xincr=tstep
      npoint=jtrflg
      asweep=atimex
      go to 40
   30 xstart=fstart
      xincr=fincr
      npoint=icalc
      asweep=afreq
c
c  construct and print the output variables with corresponding plot
c    symbols
c
   40 loct=loc+2
      if (kntr.eq.1) go to 80
      write (iofile,41)
   41 format('0legend:'/)
      do 70 i=1,kntr
      loct=loct+2
      itab(i)=nodplc(loct)
      ioutyp=nodplc(loct+1)
      itype(i)=ioutyp
      ilogy(i)=1
      if (mode.le.2) go to 50
      ilogy(i)=logopt(ioutyp)
   50 ipos=1
      call outnam(itab(i),itype(i),string,ipos)
      call move(string,ipos,ablnk,1,7)
      jstop=(ipos+6)/8
      call move(achar,1,pltsym,i,1)
      write (iofile,61) achar,(string(j),j=1,jstop)
   61 format(1x,a1,2h: ,5a8)
   70 continue
   80 if (kntr.ge.2) go to 90
      itab(1)=nodplc(loc+4)
      ioutyp=nodplc(loc+5)
      itype(1)=ioutyp
      ilogy(1)=1
      if (mode.le.2) go to 90
      ilogy(1)=logopt(ioutyp)
   90 ipos=1
      call outnam(itab(1),itype(1),string,ipos)
      call move(string,ipos,ablnk,1,7)
      jstop=(ipos+6)/8
      write (iofile,101) asweep,(string(j),j=1,jstop)
  101 format(1hx/3x,a8,4x,5a8)
      return
      end
