      subroutine runcon(id)
      implicit double precision (a-h,o-z)
c
c     this routine processes run control cards.
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
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
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
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      dimension iprnt(5),limits(4),itrlim(6),contol(6),dflts(4)
      equivalence (iprnt(1),iprnta),(limits(1),limtim),(itrlim(1),itl1),
     1   (contol(1),gmin),(dflts(1),defl)
c
c
      integer xxor
c
c  print/plot keywords
c
      dimension aopt(5)
      dimension aopts(34),lsetop(5)
      dimension aide(20)
      data aopt / 2hdc, 2htr, 2hac, 2hno, 2hdi /
c
c  options card keywords
c
      data aopts / 6hacct  , 6hlist  , 6hnomod , 6hnode  , 6hopts  ,
     1             6hitl1  , 6hitl2  , 6hitl3  , 6hitl4  , 6hitl5  ,
     2             6hitl6  , 6hlimtim, 6hlimpts, 6hlvlcod, 6hlvltim,
     3             6hgmin  , 6hreltol, 6habstol, 6hvntol , 6htrtol ,
     4             6hchgtol, 6htnom  , 6hnumdgt, 6hmaxord, 6hmethod,
     5             6hnopage, 6hmu    , 6hcptime, 6hdefl  , 6hdefw  ,
     6             6hdefad , 6hdefas , 6hpivtol, 6hpivrel /
      data lsetop / 1 ,1, 0, 1, 1 /
c
c
      data aide / 1hr,1hc,1hl,1hk,1hg,1he,1hf,1hh,1hv,1hi,1hd,1hq,1hj,
     1   1hm,1hs,1hy,1ht,4htemp,1hx,0.0d0 /
      data alsde,alsoc,alsli / 3hdec, 3hoct, 3hlin /
      data atrap, agear, auic / 4htrap, 4hgear, 3huic /
      data ablnk, ain, aout / 1h , 2hin, 3hout /
      data amiss / 8h*missing /
      data ams / 2hms /
      data minpts / 1 /
c
c
      go to (1200,1100,1650,6000,6000,1700,6000,1600,1550,2000,3600,
     1   3500,6000,1750,1300,1500,1800,4000,4100,4200,5900), id
c
c  dc transfer curves
c
 1100 ifld=2
      icvflg=0
      inum=1
 1105 anam=value(ifield+ifld)
      if(inum.gt.2) go to 6000
      id=0
      call move(anam,2,ablnk,1,7)
      if (anam.eq.aide(1)) id=1
      if (anam.eq.aide(9)) id=9
      if (anam.eq.aide(10)) id=10
      if (anam.eq.aide(17)) go to 1108
      if (id.eq.0) go to 1130
      call find(value(ifield+ifld),id,itcelm(inum),0)
      go to 1115
 1108 anam=value(ifield+ifld)
      call move(anam,5,ablnk,1,4)
      if (anam.ne.aide(18)) go to 1130
      id=18
      call find(anam,id,itcelm(inum),1)
      locs=nodplc(itcelm(inum)+1)
      nodplc(itcelm(inum)+2)=0
      value(locs)=anam
      value(locs+1)=value(itemps+1)
      call extmem(itemps,2)
      value(itemps+2)=value(itemps+1)
 1115 ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1130
      tcstar(inum)=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1130
      tcstop(inum)=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1130
      tcincr(inum)=value(ifield+ifld)
      if (tcincr(inum).eq.0.0d0) go to 1130
      temp=(tcstop(inum)-tcstar(inum))/tcincr(inum)
      if (temp.gt.0.0d0) go to 1110
      tcincr(inum)=-tcincr(inum)
      temp=-temp
 1110 itemp=idint(temp+0.5d0)+1
      itemp=max0(itemp,minpts)
      if(inum.eq.1) icvflg=itemp
      if(inum.eq.2) icvflg=itemp*icvflg
      ifld=ifld+1
      inum=2
      if(nodplc(icode+ifld)) 6000,1130,1105
 1130 write (iofile,1131)
      icvflg=0
 1131 format('0warning:  missing parameter(s) ... analysis omitted'/)
      go to 6000
c
c  frequency specification
c
 1200 ifld=2
      if (nodplc(icode+2)) 1250,1250,1210
 1210 id=0
      if (value(ifield+ifld).eq.alsde) id=1
      if (value(ifield+ifld).eq.alsoc) id=2
      if (value(ifield+ifld).eq.alsli) id=3
      if (id.eq.0) go to 1240
      idfreq=id
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1250
      if (value(ifield+ifld).le.0.0d0) go to 1250
      fincr=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1250
      if (value(ifield+ifld).le.0.0d0) go to 1250
      fstart=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1250
      if (value(ifield+ifld).le.0.0d0) go to 1250
      fstop=value(ifield+ifld)
      if (fstart.gt.fstop) go to 1260
      jacflg=fincr
      if (idfreq-2) 1215,1220,1235
 1215 fincr=dexp(xlog10/fincr)
      go to 1230
 1220 fincr=dexp(xlog2/fincr)
 1230 temp=dlog(fstop/fstart)/dlog(fincr)
      jacflg=idint(temp+0.999d0)+1
 1235 jacflg=max0(jacflg,minpts)
      if (idfreq.ne.3) go to 6000
      fincr=(fstop-fstart)/dble(max0(jacflg-1,1))
      go to 6000
 1240 write (iofile,1241) value(ifield+ifld)
 1241 format('0warning:  unknown frequency function:  ',a8,' ... analys'
     1   ,'is omitted'/)
      go to 6000
 1250 write (iofile,1251)
 1251 format('0warning:  frequency parameters incorrect ... analysis om'
     1   ,'itted'/)
      go to 6000
 1260 write (iofile,1261)
 1261 format('0warning:  start freq > stop freq ... analysis omitted'/)
      go to 6000
c
c  time specification
c
 1300 ifld=2
      if (nodplc(icode+ifld).ne.0) go to 1430
      if (value(ifield+ifld).le.0.0d0) go to 1430
      tstep=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1430
      if (value(ifield+ifld).le.0.0d0) go to 1430
      tstop=value(ifield+ifld)
      tstart=0.0d0
      delmax=tstop/50.0d0
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1310
      if (value(ifield+ifld).lt.0.0d0) go to 1430
      tstart=value(ifield+ifld)
      delmax=(tstop-tstart)/50.0d0
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1310
      if (value(ifield+ifld).le.0.0d0) go to 1430
      delmax=value(ifield+ifld)
      ifld=ifld+1
 1310 if (nodplc(icode+ifld).ne.1) go to 1320
      if (value(ifield+ifld).ne.auic) go to 1320
      nosolv=1
 1320 if (tstart.gt.tstop) go to 1440
      if (tstep.gt.tstop) go to 1430
      jtrflg=idint((tstop-tstart)/tstep+0.5d0)+1
      jtrflg=max0(jtrflg,minpts)
      go to 6000
 1430 write (iofile,1431)
 1431 format('0warning:  time parameters incorrect ... analysis omitted'
     1   /)
      go to 6000
 1440 write (iofile,1441)
 1441 format('0warning:  start time > stop time ... analysis omitted'/)
      go to 6000
c
c  transfer function
c
 1500 kssop=1
      ifld=2
      if (nodplc(icode+ifld).ne.1) go to 1530
      call outdef(ifld,1,kovar,ktype)
      if (igoof.ne.0) go to 1530
      if (ktype.ne.1) go to 1540
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 1530
      anam=value(ifield+ifld)
      call move(anam,2,ablnk,1,7)
      id=0
      if (anam.eq.aide(9)) id=9
      if (anam.eq.aide(10)) id=10
      if (id.eq.0) go to 1530
      call find(value(ifield+ifld),id,kinel,0)
      kidin=id
      go to 6000
 1530 kovar=0
      kinel=0
      write (iofile,1131)
      igoof=0
      go to 6000
 1540 kovar=0
      kinel=0
      write (iofile,1541)
 1541 format('0warning:  illegal output variable ... analysis omitted'/)
      igoof=0
      go to 6000
c
c  operating point
c
 1550 kssop=1
      go to 6000
c
c  noise analysis
c
 1600 ifld=2
      if (nodplc(icode+ifld).ne.1) go to 1610
      call outdef(ifld,2,nosout,ntype)
      if (igoof.ne.0) go to 1610
      if (ntype.ne.1) go to 1610
      if (nodplc(nosout+5).ne.0) go to 1610
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 1620
      anam=value(ifield+ifld)
      call move(anam,2,ablnk,1,7)
      id=0
      if (anam.eq.aide(9)) id=9
      if (anam.eq.aide(10)) id=10
      if (id.eq.0) go to 1620
      call find(value(ifield+ifld),id,nosin,0)
      nosprt=0
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 1605
      nosprt=dmax1(0.0d0,value(ifield+ifld))
 1605 inoise=1
      go to 6000
 1610 write (iofile,1611)
 1611 format('0warning:  voltage output unrecognizable ... analysis omit
     1ted'/)
      igoof=0
      go to 6000
 1620 write (iofile,1621)
 1621 format('0warning:  invalid input source ... analysis omitted'/)
      igoof=0
      go to 6000
c
c  distortion analysis
c
 1650 ifld=2
      if (nodplc(icode+ifld).ne.1) go to 1660
      anam=value(ifield+ifld)
      call move(anam,2,ablnk,1,7)
      if (anam.ne.aide(1)) go to 1660
      call find(value(ifield+ifld),1,idist,0)
      idprt=0
      skw2=0.9d0
      refprl=1.0d-3
      spw2=1.0d0
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 6000
      idprt=value(ifield+ifld)
      idprt=max0(idprt,0)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 6000
      if (value(ifield+ifld).le.0.001d0) go to 1670
      if (value(ifield+ifld).gt.0.999d0) go to 1670
      skw2=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 6000
      if (value(ifield+ifld).lt.1.0d-10) go to 1670
      refprl=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 6000
      if (value(ifield+ifld).lt.0.001d0) go to 1670
      spw2=value(ifield+ifld)
      go to 6000
 1660 write (iofile,1661)
 1661 format('0warning:  distortion load resistor missing ... analysis '
     1   ,'omitted'/)
      go to 6000
 1670 idist=0
      write (iofile,1671)
 1671 format('0warning:  distortion parameters incorrect ... analysis o'
     1   ,'mitted'/)
      go to 6000
c
c  fourier analysis
c
 1700 ifld=2
      if (nodplc(icode+ifld).ne.0) go to 1720
      if (value(ifield+ifld).le.0.0d0) go to 1720
      forfre=value(ifield+ifld)
 1705 ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 1710
      call outdef(ifld,2,loct,ltype)
      if (igoof.ne.0) go to 1720
      if (ltype.ne.1) go to 1720
      call extmem(ifour,1)
      nfour=nfour+1
      nodplc(ifour+nfour)=loct
      go to 1705
 1710 if (nfour.ge.1) go to 6000
 1720 write (iofile,1721)
 1721 format('0warning:  fourier parameters incorrect ... analysis omit'
     1   ,'ted'/)
      igoof=0
      nfour=0
      call clrmem(ifour)
      call getm4(ifour,0)
      go to 6000
c
c  sensitivity analysis
c
 1750 kssop=1
      ifld=1
 1760 ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 6000
      call outdef(ifld,1,loct,ltype)
      if (igoof.ne.0) go to 1780
      if (ltype.ne.1) go to 1780
      call extmem(isens,1)
      nsens=nsens+1
      nodplc(isens+nsens)=loct
      go to 1760
 1780 write (iofile,1781)
 1781 format('0warning:  output variable unrecognizable ... analysis om'
     1   ,'mitted'/)
      igoof=0
      nsens=0
      call clrmem(isens)
      call getm4(isens,0)
      go to 6000
c
c  temperature variation
c
 1800 ifld=1
 1810 ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 6000
      if (value(ifield+ifld).le.-223.0d0) go to 1810
      call extmem(itemps,1)
      numtem=numtem+1
      value(itemps+numtem)=value(ifield+ifld)
      go to 1810
c
c  options card
c
 2000 ifld=1
 2010 ifld=ifld+1
 2020 if (nodplc(icode+ifld)) 6000,2010,2030
 2030 anam=value(ifield+ifld)
      do 2040 i=1,5
      if (anam.ne.aopts(i)) go to 2040
      iprnt(i)=lsetop(i)
      ifld=ifld+1
      if(nodplc(icode+ifld).ne.0) go to 2020
      iprnt(i)=value(ifield+ifld)
      go to 2010
 2040 continue
      if (anam.eq.aopts(25)) go to 2110
      if (anam.eq.aopts(26)) go to 2120
      if (anam.eq.aopts(27)) go to 2130
      if (anam.eq.aopts(28)) go to 2150
      if (anam.eq.aopts(33)) go to 2200
      if (anam.eq.aopts(34)) go to 2250
      if (nodplc(icode+ifld+1).ne.0) go to 2510
      ifld=ifld+1
      aval=value(ifield+ifld)
      do 2050 i=6,11
      if (anam.ne.aopts(i)) go to 2050
      if(aval.le.0.0d0.and.i.ne.10) go to 2510
      itrlim(i-5)=aval
      go to 2010
 2050 continue
      if (aval.le.0.0d0) go to 2510
      do 2060 i=12,15
      if (anam.ne.aopts(i)) go to 2060
      limits(i-11)=aval
      go to 2010
 2060 continue
      do 2070 i=16,21
      if (anam.ne.aopts(i)) go to 2070
      contol(i-15)=aval
      go to 2010
 2070 continue
      do 2075 i=29,32
      if(anam.ne.aopts(i)) go to 2075
      dflts(i-28)=aval
      go to 2010
 2075 continue
      if (anam.ne.aopts(22)) go to 2080
      if (aval.lt.-223.0d0) go to 2510
      value(itemps+1)=aval
      go to 2010
 2080 if (anam.ne.aopts(23)) go to 2100
      ndigit=aval
      if (ndigit.le.7) go to 2090
      ndigit=7
      write (iofile,2081) ndigit
 2081 format('0warning:  numdgt may not exceed',i2,
     1 ';  maximum value assumed'/)
 2090 numdgt=ndigit
      go to 2010
 2100 if (anam.ne.aopts(24)) go to 2500
      n=aval
      if ((n.le.1).or.(n.ge.7)) go to 2510
      maxord=n
      go to 2010
 2110 if (nodplc(icode+ifld+1).ne.1) go to 2510
      ifld=ifld+1
      anam=value(ifield+ifld)
      call move(anam,5,ablnk,1,4)
      jtype=0
      if (anam.eq.atrap) jtype=1
      if (anam.eq.agear) jtype=2
      if (jtype.eq.0) go to 2510
      method=jtype
      go to 2010
 2120 nopage=1
      go to 2010
 2130 ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,2140,2030
 2140 aval=value(ifield+ifld)
      if(aval.lt.0.0d0.or.aval.gt.0.500001d0) go to 2510
      xmu=aval
      go to 2010
 2150 ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,2160,2030
 2160 aval=value(ifield+ifld)
      maxtim=aval
      go to 2010
 2200 ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,2210,2030
 2210 aval=value(ifield+ifld)
      if (aval.gt.1.0d0) go to 2510
      pivtol=aval
      go to 2010
 2250 ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,2260,2030
 2260 aval=value(ifield+ifld)
      if (aval.gt.1.0d0) go to 2510
      pivrel=aval
      go to 2010
 2500 write (iofile,2501) anam
 2501 format('0warning:  unknown option:  ',a8,' ... ignored'/)
      go to 2010
 2510 write (iofile,2511) anam
 2511 format('0warning:  illegal value specified for option:  ',a8,' ...
     1 ignored'/)
      go to 2010
c
c  print card
c
 3500 iprpl=0
      go to 3610
c
c  plot (and print) card
c
 3600 iprpl=1
 3610 ifld=2
 3613 anam=amiss
      if (nodplc(icode+ifld).ne.1) go to 3950
      anam=value(ifield+ifld)
      ms=0
      if (xxor(anam,ams).ne.0) go to 3615
      ms=1
      ifld=3
      if (nodplc(icode+ifld).ne.1) go to 3970
      anam=value(ifield+ifld)
 3615 call move(anam,3,ablnk,1,6)
      do 3620 i=1,5
      if (anam.ne.aopt(i)) go to 3620
      ktype=i
      go to 3630
 3620 continue
      go to 3950
 3630 id=30+5*iprpl+ktype
      call find(dble(jelcnt(id)),id,loc,1)
      nodplc(loc+2)=ktype
      if (ms.eq.0) go to 3635
      locv=nodplc(loc+1)
      value(locv)=0.0d0
 3635 numout=0
 3640 ifld=ifld+1
      if (nodplc(icode+ifld)) 3900,3640,3650
 3650 call outdef(ifld,ktype,loct,ltype)
      if (igoof.ne.0) go to 3970
      if (iprpl.eq.0) go to 3660
      plimlo=0.0d0
      plimhi=0.0d0
      if (nodplc(icode+ifld+1).ne.0) go to 3660
      if (nodplc(icode+ifld+2).ne.0) go to 3660
      plimlo=value(ifield+ifld+1)
      plimhi=value(ifield+ifld+2)
      ifld=ifld+2
 3660 numout=numout+1
      lspot=loc+2*numout+2
      nodplc(lspot)=loct
      nodplc(lspot+1)=ltype
      if (iprpl.eq.0) go to 3670
      locv=nodplc(loc+1)
      lspot=locv+2*numout-1
      value(lspot)=plimlo
      value(lspot+1)=plimhi
 3670 if (numout.eq.8) go to 3900
      go to 3640
 3900 nodplc(loc+3)=numout
      if (iprpl.eq.0) go to 6000
c...  propogate plot limits downward
      if (numout.le.1) go to 6000
      locv=nodplc(loc+1)
      lspot=locv+2*numout-1
      plimlo=value(lspot)
      plimhi=value(lspot+1)
      i=numout-1
 3905 lspot=lspot-2
      if (value(lspot).ne.0.0d0) go to 3910
      if (value(lspot+1).ne.0.0d0) go to 3910
      value(lspot)=plimlo
      value(lspot+1)=plimhi
      go to 3920
 3910 plimlo=value(lspot)
      plimhi=value(lspot+1)
 3920 i=i-1
      if (i.ge.1) go to 3905
      go to 6000
c
c     errors
c
 3950 write (iofile,3951) anam
 3951 format('0warning:  unknown analysis mode:  ',a8,
     1  ' ... line ignored'/)
      go to 6000
 3970 write (iofile,3971)
 3971 format('0warning:  unrecognizable output variable on above line'/)
      igoof=0
      go to 3640
c
c  width card
c
 4000 ifld=1
 4010 ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 6000
 4020 anam=value(ifield+ifld)
      if (anam.ne.ain) go to 4040
      ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,4030,4020
 4030 iwidth=value(ifield+ifld)
      iwidth=min0(max0(iwidth,10),120)
      go to 4010
 4040 if (anam.ne.aout) go to 6000
      ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,4050,4020
 4050 lwidth=dmin1(dmax1(value(ifield+ifld),72.0d0),132.0d0)
      go to 4010
c
c  nodeset statement
c
 4100 ifld=1
 4110 ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,4120,4110
 4120 nodnum=value(ifield+ifld)
      if(nodnum.le.0) go to 4190
      ifld=ifld+1
      if(nodplc(icode+ifld)) 4180,4130,4170
 4130 call sizmem(nsnod,nic)
      call extmem(nsnod,1)
      call extmem(nsval,1)
      nodplc(nsnod+nic+1)=nodnum
      value(nsval+nic+1)=value(ifield+ifld)
      go to 4110
c
c  errors on .nodeset statement
c
 4170 write(iofile,4171) value(ifield+ifld)
 4171 format('0warning: out-of-place non-numeric field ',a8,
     1 ' skipped'/)
      go to 4110
 4180 write(iofile,4181) nodnum
 4181 format('0warning: initial value missing for node ',i5,/)
      go to 6000
 4190 write(iofile,4191)
 4191 format('0warning: attempt to specify initial condition for ',
     1 'ground ingnored',/)
      ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,4110,4170
c
c  initial conditions statement
c
 4200 ifld=1
 4210 ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,4220,4210
 4220 nodnum=value(ifield+ifld)
      if(nodnum.le.0) go to 4290
      ifld=ifld+1
      if(nodplc(icode+ifld)) 4280,4230,4270
 4230 call sizmem(icnod,nic)
      call extmem(icnod,1)
      call extmem(icval,1)
      nodplc(icnod+nic+1)=nodnum
      value(icval+nic+1)=value(ifield+ifld)
      go to 4210
c
c  errors on .ic statement
c
 4270 write(iofile,4271) value(ifield+ifld)
 4271 format('0warning: out-of-place non-numeric field ',a8,
     1 ' skipped'/)
      go to 4210
 4280 write(iofile,4281) nodnum
 4281 format('0warning: initial value missing for node ',i5,/)
      go to 6000
 4290 write(iofile,4291)
 4291 format('0warning: attempt to specify initial condition for ',
     1 'ground ignored',/)
      ifld=ifld+1
      if(nodplc(icode+ifld)) 6000,4210,4270
c
c     :debug: statement
c     sample debug line: .:debug: 5=3 17=5
c
 5900 ifld=1
 5910 ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,5920,5910
 5920 index=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld)) 6000,5930,5910
 5930 ival=value(ifield+ifld)
      if (index.lt.1) go to 5910
      if (index.gt.20) go to 5910
      write(iofile,5931) index,ival
 5931 format(' *debug*:  runcon - idebug(',i2,') set to ',i10)
      idebug(index)=ival
      go to 5910
c
c  finished
c
 6000 return
      end
