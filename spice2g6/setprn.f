      subroutine setprn(loc)
      implicit double precision (a-h,o-z)
c
c     this routine formats the column headers for tabular listings of
c output variables.
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
      data ablnk, atimex, afreq / 1h , 6h  time, 6h  freq /
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
c  construct and print the output variable names
c
   40 loct=loc+2
      ipos=1
      npos=ipos+numdgt+8
      do 90 i=1,kntr
      loct=loct+2
      itab(i)=nodplc(loct)
      itype(i)=nodplc(loct+1)
      call outnam(itab(i),itype(i),string,ipos)
      if (ipos.ge.npos) go to 70
      do 60 j=ipos,npos
      call move(string,j,ablnk,1,1)
   60 continue
      ipos=npos
      go to 80
   70 call move(string,ipos,ablnk,1,1)
      ipos=ipos+1
   80 npos=npos+numdgt+8
   90 continue
      call move(string,ipos,ablnk,1,7)
      jstop=(ipos+6)/8
      write (iofile,91) asweep,(string(j),j=1,jstop)
   91 format(/3x,a8,5x,14a8,a4)
      write (iofile,101)
  101 format(1hx/1h )
      return
      end
