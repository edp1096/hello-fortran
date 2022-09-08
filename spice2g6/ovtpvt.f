c spice version 2g.6  sccsid=ovtpvt.ma 3/15/83
      subroutine ovtpvt
      implicit double precision (a-h,o-z)
c
c
c     this routine generates the requested tabular listings of analysis
c results.  it calls plot to generate line-printer plots.
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
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
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
      complex cval
      dimension prform(3)
      dimension subtit(4,3)
      data subtit / 8hdc trans, 8hfer curv, 8hes      , 8h        ,
     1              8htransien, 8ht analys, 8his      , 8h        ,
     2              8hac analy, 8hsis     , 8h        , 8h         /
      data prform / 8h(1pe11.3, 8h,2x,8e00, 8h.00)     /
      data aper,rprn / 1h., 1h) /
c
      call second(t1)
      if (icalc.le.0) go to 1000
      call crunch
      if (nogo.lt.0) go to 1000
c
c  construct format statement to be used for printing the outputs
c
      ifract=max0(numdgt-1,0)
      ifwdth=ifract+9
      ipos=15
      call alfnum(ifwdth,prform,ipos)
      call move(prform,ipos,aper,1,1)
      ipos=ipos+1
      call alfnum(ifract,prform,ipos)
      call move(prform,ipos,rprn,1,1)
c
      noprln=min0(8,(lwidth-12)/ifwdth)
      if (mode-2) 50,60,300
   50 numout=jelcnt(41)+1
      go to 70
   60 numout=jelcnt(42)+1
c
c  dc and transient analysis printing
c
   70 loc=locate(30+mode)
   80 if (loc.eq.0) go to 200
      kntr=min0(noprln,nodplc(loc+3))
      if (kntr.le.0) go to 120
      call title(1,lwidth,1,subtit(1,mode))
      call setprn(loc)
c
c  get buffer space
c
      call getm8(locx,npoint)
      call getm8(locy,kntr*npoint)
c
c  interpolate outputs
c
      call ntrpl8(locx,locy,numpnt)
c
c  print outputs
c
      do 100 i=1,numpnt
      xvar=value(locx+i)
      locyt=locy
      do 90 k=1,kntr
      yvar(k)=value(locyt+i)
      locyt=locyt+npoint
   90 continue
      write (iofile,prform) xvar,(yvar(k),k=1,kntr)
  100 continue
      write (iofile,111)
  111 format(1hy)
      call clrmem(locx)
      call clrmem(locy)
  120 loc=nodplc(loc)
      go to 80
c
c  dc and transient analysis plotting
c
  200 loc=locate(35+mode)
  210 if (loc.eq.0) go to 250
      kntr=nodplc(loc+3)
      if (kntr.le.0) go to 220
      locv=nodplc(loc+1)
      call title(1,lwidth,1,subtit(1,mode))
      call setplt(loc)
c
c     get buffer space
c
      call getm8(locx,npoint)
      call getm8(locy,kntr*npoint)
c
c  interpolate outputs and load plot buffers
c
      call ntrpl8(locx,locy,numpnt)
      call plot(numpnt,locx,locy,locv)
      call clrmem(locx)
      call clrmem(locy)
  220 loc=nodplc(loc)
      go to 210
c
c  fourier analysis
c
  250 if (mode.eq.1) go to 1000
      if (nfour.eq.0) go to 1000
      if (nogo.ne.0) go to 1000
      call fouran
      go to 1000
c
c  ac analysis printing
c
  300 numout=jelcnt(43)+jelcnt(44)+jelcnt(45)+1
      do 599 id=33,35
      loc=locate(id)
  320 if (loc.eq.0) go to 599
      kntr=min0(noprln,nodplc(loc+3))
      if (kntr.le.0) go to 595
      call title(1,lwidth,1,subtit(1,mode))
      call setprn(loc)
c
c  print ac outputs
c
      lout=loutpt
      do 590 i=1,icalc
      xvar=dble(real(cvalue(lout+1)))
      do 500 k=1,kntr
      iseq=itab(k)
      iseq=nodplc(iseq+4)
      cval=cvalue(lout+iseq)
      ktype=itype(k)
      go to (450,450,430,440,450,450), ktype
  430 yvar(k)=dble(real(cval))
      go to 500
  440 yvar(k)=dble(aimag(cval))
      go to 500
  450 call magphs(cval,xmag,xphs)
      go to (460,460,430,440,470,465), ktype
  460 yvar(k)=xmag
      go to 500
  465 yvar(k)=20.0d0*dlog10(xmag)
      go to 500
  470 yvar(k)=xphs
  500 continue
      lout=lout+numout
  580 write (iofile,prform) xvar,(yvar(k),k=1,kntr)
  590 continue
      write (iofile,111)
  595 loc=nodplc(loc)
      go to 320
  599 continue
c
c  ac analysis plotting
c
      do 760 id=38,40
      loc=locate(id)
  610 if (loc.eq.0) go to 760
      kntr=nodplc(loc+3)
      if (kntr.le.0) go to 750
      locv=nodplc(loc+1)
      call title(1,lwidth,1,subtit(1,mode))
      call setplt(loc)
c
      call getm8(locx,icalc)
      call getm8(locy,kntr*icalc)
c
c     load plot buffers
c
      lout=loutpt
      do 710 i=1,icalc
      xvar=dble(real(cvalue(lout+1)))
      locyt=locy
      do 700 k=1,kntr
      iseq=itab(k)
      iseq=nodplc(iseq+4)
      cval=cvalue(lout+iseq)
      ktype=itype(k)
      go to (670,670,650,660,670,670), ktype
  650 yvr=dble(real(cval))
      go to 695
  660 yvr=dble(aimag(cval))
      go to 695
  670 call magphs(cval,xmag,xphs)
      go to (680,680,650,660,690,685), ktype
  680 yvr=dlog10(xmag)
      go to 695
  685 yvr=20.0d0*dlog10(xmag)
      go to 695
  690 yvr=xphs
  695 value(locyt+i)=yvr
      locyt=locyt+icalc
  700 continue
      value(locx+i)=xvar
      lout=lout+numout
  710 continue
      call plot(icalc,locx,locy,locv)
      call clrmem(locx)
      call clrmem(locy)
  750 loc=nodplc(loc)
      go to 610
  760 continue
c
c  finished
c
 1000 call clrmem(loutpt)
      call second(t2)
      rstats(11)=rstats(11)+t2-t1
      return
      end
