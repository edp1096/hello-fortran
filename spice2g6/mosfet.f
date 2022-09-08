      subroutine mosfet
      implicit double precision (a-h,o-z)
c
c     this routine processes mosfets for dc and transient analyses.
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
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension vbdo(1),vbso(1),vgso(1),vdso(1),cdo(1),cbso(1),cbdo(1),
     1   gmo(1),gdso(1),gmbso(1),gbdo(1),gbso(1),
     2   qb(1),cqb(1),qg(1),cqg(1),qd(1),cqd(1),
     3   cggbo(1),cgdbo(1),cgsbo(1),cbgbo(1),cbdbo(1),cbsbo(1),
     4   cgbo(1),cgdo(1),cgso(1),vono(1),vdsato(1)
      dimension qbd(1),cqbd(1),qbs(1),cqbs(1),
     1          qgs(1),cqgs(1),qgd(1),cqgd(1),qgb(1),cqgb(1)
      equivalence (vbdo (1),value( 1)),(vbso (1),value( 2)),
     1            (vgso (1),value( 3)),(vdso (1),value( 4)),
     2            (cdo  (1),value( 5)),(cbso (1),value( 6)),
     3            (cbdo (1),value( 7)),(gmo  (1),value( 8)),
     4            (gdso (1),value( 9)),(gmbso(1),value(10)),
     5            (gbdo (1),value(11)),(gbso (1),value(12)),
     6            (qb   (1),qgs  ( 1), value(13)),
     7            (cqb  (1),cqgs ( 1), value(14)),
     8            (qg   (1),qgd  ( 1), value(15)),
     9            (cqg  (1),cqgd ( 1), value(16)),
     a            (qd   (1),qgb  ( 1), value(17)),
     b            (cqd  (1),cqgb ( 1), value(18)),
     c            (cggbo(1),cgbo  (1), value(19)),
     d            (cgdbo(1),cgdo  (1), value(20)),
     e            (cgsbo(1),cgso  (1), value(21)),
     f            (cbgbo(1),vono  (1), value(22)),
     g            (cbdbo(1),vdsato(1), value(23)),
     h            (cbsbo(1),           value(24))
      equivalence (qbd  (1),value(25)),(cqbd (1),value(26)),
     1            (qbs  (1),value(27)),(cqbs (1),value(28))
c
c
      loc=locate(14)
   10 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) return
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      locm=nodplc(loc+8)
      ioff=nodplc(loc+9)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=nodplc(loc+26)
c
c  dc model parameters
c
      xj=value(locm+27)
      xld=value(locm+28)
      xl=value(locv+1)-2.0d0*xld
      xw=value(locv+2)
      devmod=value(locv+8)
      vto=type*value(locm+2)
      vdsat=type*value(locv+10)
      vinit=value(locm+43)
      ad=value(locv+3)
      as=value(locv+4)
      pd=value(locv+11)
      ps=value(locv+12)
      if (value(locm+21).eq.0.0d0.
     1   or.ad.eq.0.0d0.or.as.eq.0.0d0) go to 12
      cdsat=value(locm+21)*ad
      cssat=value(locm+21)*as
      go to 15
   12 cdsat=value(locm+11)
      cssat=value(locm+11)
   15 if ((value(locm+7).le.0.0d0).and.
     1              (value(locm+8).le.0.0d0)) go to 17
      gdpr=value(locm+7)
      gspr=value(locm+8)
      go to 19
   17 gdpr=value(locm+16)/value(locv+13)
      gspr=value(locm+16)/value(locv+14)
   19 covlgs=value(locm+13)*xw
      covlgd=value(locm+14)*xw
      covlgb=value(locm+15)*xl
      lev=value(locm+1)
c
c     mos1, mos2 and mos3 model parameters
c
      beta=value(locm+3)*xw/xl
      gamma=value(locm+4)
      phi=value(locm+5)
      xlamda=value(locm+6)
      phib=value(locm+12)
      cox=value(locm+22)*xw*xl
      xnsub=value(locm+23)
      xnfs=value(locm+25)
      uo=value(locm+29)
      vbp=value(locm+30)
      uexp=value(locm+31)
      utra=value(locm+32)
      vbi=type*value(locm+44)
      xd=value(locm+45)
      vmax=value(locm+33)
      xneff=value(locm+34)
      xqco=value(locm+35)
      fnarrw=value(locm+39)
      if (lev.eq.3) fnarrw=fnarrw/xw
c
c     initialization
c
      icheck=1
      ibypas=0
      go to (100,20,30,50,60,70), initf
   20 if (ioff.ne.0) go to 40
      vds=type*value(locv+5)
      vgs=type*value(locv+6)
      vbs=type*value(locv+7)
      if (vds.ne.0.0d0) go to 300
      if (vgs.ne.0.0d0) go to 300
      if (vbs.ne.0.0d0) go to 300
      if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 300
      vbs=vinit
      vgs=vto
      vds=0.0d0
      go to 300
   30 if (ioff.eq.0) go to 100
   40 vbs=0.0d0
      vgs=0.0d0
      vds=0.0d0
      go to 300
   50 vbs=vbso(lx0+loct)
      vgs=vgso(lx0+loct)
      vds=vdso(lx0+loct)
      go to 300
   60 vbs=vbso(lx1+loct)
      vgs=vgso(lx1+loct)
      vds=vdso(lx1+loct)
      go to 300
   70 xfact=delta/delold(2)
      vbso(lx0+loct)=vbso(lx1+loct)
      vbs=(1.0d0+xfact)*vbso(lx1+loct)-xfact*vbso(lx2+loct)
      vgso(lx0+loct)=vgso(lx1+loct)
      vgs=(1.0d0+xfact)*vgso(lx1+loct)-xfact*vgso(lx2+loct)
      vdso(lx0+loct)=vdso(lx1+loct)
      vds=(1.0d0+xfact)*vdso(lx1+loct)-xfact*vdso(lx2+loct)
      vbdo(lx0+loct)=vbso(lx0+loct)-vdso(lx0+loct)
      cdo(lx0+loct)=cdo(lx1+loct)
      cbso(lx0+loct)=cbso(lx1+loct)
      cbdo(lx0+loct)=cbdo(lx1+loct)
      gmo(lx0+loct)=gmo(lx1+loct)
      gdso(lx0+loct)=gdso(lx1+loct)
      gmbso(lx0+loct)=gmbso(lx1+loct)
      gbdo(lx0+loct)=gbdo(lx1+loct)
      gbso(lx0+loct)=gbso(lx1+loct)
      cggbo(lx0+loct)=cggbo(lx1+loct)
      cgdbo(lx0+loct)=cgdbo(lx1+loct)
      cgsbo(lx0+loct)=cgsbo(lx1+loct)
      cbgbo(lx0+loct)=cbgbo(lx1+loct)
      cbdbo(lx0+loct)=cbdbo(lx1+loct)
      cbsbo(lx0+loct)=cbsbo(lx1+loct)
      go to 110
c
c  compute new nonlinear branch voltages
c
  100 vbs=type*(value(lvnim1+node4)-value(lvnim1+node6))
      vgs=type*(value(lvnim1+node2)-value(lvnim1+node6))
      vds=type*(value(lvnim1+node5)-value(lvnim1+node6))
  110 vbd=vbs-vds
      vgd=vgs-vds
      vgdo=vgso(lx0+loct)-vdso(lx0+loct)
      delvbs=vbs-vbso(lx0+loct)
      delvbd=vbd-vbdo(lx0+loct)
      delvgs=vgs-vgso(lx0+loct)
      delvds=vds-vdso(lx0+loct)
      delvgd=vgd-vgdo
      if (devmod.lt.0.0d0) go to 120
      cdhat=cdo(lx0+loct)-gbdo(lx0+loct)*delvbd+gmbso(lx0+loct)*delvbs
     1   +gmo(lx0+loct)*delvgs+gdso(lx0+loct)*delvds
      go to 130
  120 cdhat=cdo(lx0+loct)-(gbdo(lx0+loct)-gmbso(lx0+loct))*delvbd
     1   -gmo(lx0+loct)*delvgd+gdso(lx0+loct)*delvds
  130 cbhat=cbso(lx0+loct)+cbdo(lx0+loct)+gbdo(lx0+loct)*delvbd
     1   +gbso(lx0+loct)*delvbs
c
c  bypass if solution has not changed
c
      if (initf.eq.6) go to 200
      tol=reltol*dmax1(dabs(vbs),dabs(vbso(lx0+loct)))+vntol
      if (dabs(delvbs).ge.tol) go to 200
      tol=reltol*dmax1(dabs(vbd),dabs(vbdo(lx0+loct)))+vntol
      if (dabs(delvbd).ge.tol) go to 200
      tol=reltol*dmax1(dabs(vgs),dabs(vgso(lx0+loct)))+vntol
      if (dabs(delvgs).ge.tol) go to 200
      tol=reltol*dmax1(dabs(vds),dabs(vdso(lx0+loct)))+vntol
      if (dabs(delvds).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cdhat),dabs(cdo(lx0+loct)))+abstol
      if (dabs(cdhat-cdo(lx0+loct)).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cbhat),dabs(cbso(lx0+loct)+cbdo(lx0+loct)))
     1   +abstol
      if (dabs(cbhat-(cbso(lx0+loct)+cbdo(lx0+loct))).ge.tol) go to 200
      vbs=vbso(lx0+loct)
      vbd=vbdo(lx0+loct)
      vgs=vgso(lx0+loct)
      vds=vdso(lx0+loct)
      vgd=vgs-vds
      vgb=vgs-vbs
      cd=cdo(lx0+loct)
      cbs=cbso(lx0+loct)
      cbd=cbdo(lx0+loct)
      cdrain=devmod*(cd+cbd)
      gm=gmo(lx0+loct)
      gds=gdso(lx0+loct)
      gmbs=gmbso(lx0+loct)
      gbd=gbdo(lx0+loct)
      gbs=gbso(lx0+loct)
      devmod=value(locv+8)
      if (mode.ne.1) go to 135
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 135
      if (xqco.gt.0.5d0) go to 742
      go to 850
  135 if (xqco.le.0.5d0) go to 140
      cgb=cgbo(lx0+loct)
      cgd=cgdo(lx0+loct)
      cgs=cgso(lx0+loct)
      vgs1=vgso(lx1+loct)
      vgb1=vgs1-vbso(lx1+loct)
      vgd1=vgs1-vdso(lx1+loct)
      go to 735
  140 cggb=cggbo(lx0+loct)
      cgdb=cgdbo(lx0+loct)
      cgsb=cgsbo(lx0+loct)
      cbgb=cbgbo(lx0+loct)
      cbdb=cbdbo(lx0+loct)
      cbsb=cbsbo(lx0+loct)
      xqc=value(locv+15)
      ibypas=1
      go to 755
c
c  limit nonlinear branch voltages
c
  200 von=type*value(locv+9)
      if (vdso(lx0+loct).lt.0.0d0) go to 205
      call fetlim(vgs,vgso(lx0+loct),von)
      vds=vgs-vgd
      call limvds(vds,vdso(lx0+loct))
      vgd=vgs-vds
      go to 210
  205 call fetlim(vgd,vgdo,von)
      vds=vgs-vgd
      call limvds(-vds,-vdso(lx0+loct))
      vgs=vgd+vds
  210 if (vds.lt.0.0d0) go to 220
      vcrit=vt*dlog(vt/(root2*cssat))
      call pnjlim(vbs,vbso(lx0+loct),vt,vcrit,icheck)
      vbd=vbs-vds
      go to 300
  220 vcrit=vt*dlog(vt/(root2*cdsat))
      call pnjlim(vbd,vbdo(lx0+loct),vt,vcrit,icheck)
      vbs=vbd+vds
c
c  determine dc current and derivatives
c
  300 vbd=vbs-vds
      vgd=vgs-vds
      vgb=vgs-vbs
      if (vbs.gt.0.0d0) go to 310
      gbs=cssat/vt
      cbs=gbs*vbs
      gbs=gbs+gmin
      go to 320
  310 evbs=dexp(vbs/vt)
      gbs=cssat*evbs/vt+gmin
      cbs=cssat*(evbs-1.0d0)
  320 if (vbd.gt.0.0d0) go to 330
      gbd=cdsat/vt
      cbd=gbd*vbd
      gbd=gbd+gmin
      go to 400
  330 evbd=dexp(vbd/vt)
      gbd=cdsat*evbd/vt+gmin
      cbd=cdsat*(evbd-1.0d0)
c
c  compute drain current and derivatives
c
  400 if (vds.lt.0.0d0) go to 450
c
c  normal mode
c
      devmod=1.0d0
      value(locv+8)=devmod
      go to (405,410,415), lev
  405 call moseq1(vds,vbs,vgs,gm,gds,gmbs)
      go to 460
  410 call moseq2(vds,vbs,vgs,gm,gds,gmbs,
     1   qgate,qchan,qbulk,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      go to 460
  415 call moseq3(vds,vbs,vgs,gm,gds,gmbs,
     1   qgate,qchan,qbulk,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      go to 460
c
c  inverse mode
c
  450 devmod=-1.0d0
      value(locv+8)=devmod
      go to (452,453,454), lev
  452 call moseq1(-vds,vbd,vgd,gm,gds,gmbs)
      go to 460
  453 call moseq2(-vds,vbd,vgd,gm,gds,gmbs,
     1   qgate,qchan,qbulk,cggb,cgsb,cgdb,cbgb,cbsb,cbdb)
      go to 460
  454 call moseq3(-vds,vbd,vgd,gm,gds,gmbs,
     1   qgate,qchan,qbulk,cggb,cgsb,cgdb,cbgb,cbsb,cbdb)
  460 value(locv+9)=type*von
      value(locv+10)=type*vdsat
      if (xqco.le.0.5d0) value(locv+15)=xqc
c
c  compute equivalent drain current source
c
  490 cd=devmod*cdrain-cbd
      if (mode.ne.1) go to 500
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 500
      if (initf.eq.4) go to 500
      go to 650
c
c  charge storage elements
c
c.. bulk-drain and bulk-source depletion capacitances
c
  500 czbd=0.0d0
      czbs=0.0d0
      czbdsw=0.0d0
      czbssw=0.0d0
      if ((value(locm+9).eq.0.0d0).or.(value(locm+10).eq.0.0d0))
     1   go to 505
      czbd=value(locm+9)
      czbs=value(locm+10)
      go to 510
  505 if (value(locm+17).eq.0.0d0) go to 510
      czbd=value(locm+17)*ad
      czbs=value(locm+17)*as
  510 if (value(locm+19).eq.0.0d0) go to 515
      czbdsw=value(locm+19)*pd
      czbssw=value(locm+19)*ps
  515 phib=value(locm+12)
      xmj=value(locm+18)
      xmjsw=value(locm+20)
      twop=phib+phib
      fcpb=value(locm+38)
      fcpb2=fcpb*fcpb
      f1=value(locm+40)
      f2=value(locm+41)
      f3=value(locm+42)
      czsf2=czbs/f2
      czswf2=czbssw/f2
      czdf2=czbd/f2
      czdwf2=czbdsw/f2
      if (vbs.ge.fcpb) go to 520
      arg=1.0d0-vbs/phib
      sarg=dexp(-xmj*dlog(arg))
      sargsw=dexp(-xmjsw*dlog(arg))
      qbs(lx0+loct)=phib*(czbs*(1.0d0-arg*sarg)/(1.0d0-xmj)
     1                +czbssw*(1.0d0-arg*sargsw)/(1.0d0-xmjsw))
      capbs=czbs*sarg+czbssw*sargsw
      go to 525
  520 qbs(lx0+loct)=f1*(czbs+czbssw)+f3*(vbs-fcpb)*(czsf2+czswf2)
     1    +(vbs*vbs-fcpb*fcpb)*(czsf2*xmj+czswf2*xmjsw)
      capbs=f3*(czsf2+czswf2)+vbs/phib*(czsf2*xmj+czswf2*xmjsw)
  525 if (vbd.ge.fcpb) go to 530
      arg=1.0d0-vbd/phib
      sarg=dexp(-xmj*dlog(arg))
      sargsw=dexp(-xmjsw*dlog(arg))
      qbd(lx0+loct)=phib*(czbd*(1.0d0-arg*sarg)/(1.0d0-xmj)
     1              +czbdsw*(1.0d0-arg*sargsw)/(1.0d0-xmjsw))
      capbd=czbd*sarg+czbdsw*sargsw
      go to 560
  530 qbd(lx0+loct)=f1*(czbd+czbdsw)+f3*(vbd-fcpb)*(czdf2+czdwf2)
     1    +(vbd*vbd-fcpb*fcpb)*(czdf2*xmj+czdwf2*xmjsw)
      capbd=f3*(czdf2+czdwf2)+vbd/phib*(czdf2*xmj+czdwf2*xmjsw)
c
  560 if (xqco.le.0.5d0) go to 650
      if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 650
      if (initf.ne.4) go to 600
      go to 705
c
cc    calculate equivalent conductances and currents for
cc    depletion capacitors
c
  600 if (initf.ne.5) go to 610
      qbd(lx1+loct)=qbd(lx0+loct)
      qbs(lx1+loct)=qbs(lx0+loct)
  610 call intgr8(geq,ceq,capbd,loct+24)
      gbd=gbd+geq
      cbd=cbd+cqbd(lx0+loct)
      cd=cd-cqbd(lx0+loct)
      call intgr8(geq,ceq,capbs,loct+26)
      gbs=gbs+geq
      cbs=cbs+cqbs(lx0+loct)
      if (initf.ne.5) go to 650
      cqbd(lx1+loct)=cqbd(lx0+loct)
      cqbs(lx1+loct)=cqbs(lx0+loct)
c
c  check convergence
c
  650 if (initf.ne.3) go to 660
      if (ioff.ne.0) go to 680
  660 if (icheck.eq.1) go to 670
      tol=reltol*dmax1(dabs(cdhat),dabs(cd))+abstol
      if (dabs(cdhat-cd).ge.tol) go to 670
      tol=reltol*dmax1(dabs(cbhat),dabs(cbs+cbd))+abstol
      if (dabs(cbhat-(cbs+cbd)).le.tol) go to 680
  670 noncon=noncon+1
  680 vbso(lx0+loct)=vbs
      vbdo(lx0+loct)=vbd
      vgso(lx0+loct)=vgs
      vdso(lx0+loct)=vds
      cdo(lx0+loct)=cd
      cbso(lx0+loct)=cbs
      cbdo(lx0+loct)=cbd
      gmo(lx0+loct)=gm
      gdso(lx0+loct)=gds
      gmbso(lx0+loct)=gmbs
      gbdo(lx0+loct)=gbd
      gbso(lx0+loct)=gbs
      if (xqco.le.0.5d0) go to 690
      vono(lx0+loct)=von
      vdsato(lx0+loct)=vdsat
      go to 700
  690 cggbo(lx0+loct)=cggb
      cgdbo(lx0+loct)=cgdb
      cgsbo(lx0+loct)=cgsb
      cbgbo(lx0+loct)=cbgb
      cbdbo(lx0+loct)=cbdb
      cbsbo(lx0+loct)=cbsb
      go to 750
c
c     xqco > 0.5d0 use meyer"s capacitor model
c
  700 if (mode.ne.1) go to 705
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 705
      if (initf.eq.4) go to 705
      go to 742
c
c     calculate meyer's capacitors
c
  705 von1=von
      vgs1=vgs
      vgd1=vgd
      vgb1=vgs-vbs
      vdsat1=vdsat
      if ((mode.ne.2).or.(initf.eq.5)) go to 710
      von1=vono(lx1+loct)
      vgs1=vgso(lx1+loct)
      vgd1=vgs1-vdso(lx1+loct)
      vgb1=vgs1-vbso(lx1+loct)
      vdsat1=vdsato(lx1+loct)
  710 if (devmod.lt.0.0d0) go to 715
      call cmeyer (vgs1,vgd1,vgb1,von1,vdsat1,vgs,vgd,vgb,
     1   covlgs,covlgd,covlgb,cgs1,cgd1,cgb1,cgs,cgd,cgb)
      go to 720
  715 call cmeyer (vgd1,vgs1,vgb1,von1,vdsat1,vgd,vgs,vgb,
     1   covlgd,covlgs,covlgb,cgd1,cgs1,cgb1,cgd,cgs,cgb)
  720 cgs=0.5d0*(cgs+cgs1)
      cgd=0.5d0*(cgd+cgd1)
      cgb=0.5d0*(cgb+cgb1)
c
c     store small-signal parameters (for meyer"s model)
c
      if (mode.ne.1) go to 730
      if (initf.ne.4) go to 730
      value(lx0+loct+24)=capbd
      value(lx0+loct+26)=capbs
      value(lx0+loct+12)=cgs-covlgs
      value(lx0+loct+14)=cgd-covlgd
      value(lx0+loct+16)=cgb-covlgb
      go to 1000
cc
  730 if (initf.ne.6) go to 735
      qgs(lx0+loct)=(1.0d0+xfact)*qgs(lx1+loct)-xfact*qgs(lx2+loct)
      qgd(lx0+loct)=(1.0d0+xfact)*qgd(lx1+loct)-xfact*qgd(lx2+loct)
      qgb(lx0+loct)=(1.0d0+xfact)*qgb(lx1+loct)-xfact*qgb(lx2+loct)
      go to 745
  735 qgs(lx0+loct)=(vgs-vgs1)*cgs
      qgd(lx0+loct)=(vgd-vgd1)*cgd
      qgb(lx0+loct)=(vgb-vgb1)*cgb
      if((mode.ne.2).or.(initf.eq.5))go to 740
      qgs(lx0+loct)=qgs(lx0+loct)+qgs(lx1+loct)
      qgd(lx0+loct)=qgd(lx0+loct)+qgd(lx1+loct)
      qgb(lx0+loct)=qgb(lx0+loct)+qgb(lx1+loct)
  740 if((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 742
      if (initf.ne.5) go to 745
      qgs(lx0+loct)=cgs*vgs
      qgd(lx0+loct)=cgd*vgd
      qgb(lx0+loct)=cgb*vgb
      qgs(lx1+loct)=qgs(lx0+loct)
      qgd(lx1+loct)=qgd(lx0+loct)
      qgb(lx1+loct)=qgb(lx0+loct)
c
c     initialize to zero charge conductances and current
c
  742 gcgs=0.0d0
      ceqgs=0.0d0
      gcgd=0.0d0
      ceqgd=0.0d0
      gcgb=0.0d0
      ceqgb=0.0d0
      go to  870
cc
 745  if(cgs.eq.0.0d0) value(lx0+loct+13)=0.0d0
      if(cgd.eq.0.0d0) value(lx0+loct+15)=0.0d0
      if(cgb.eq.0.0d0) value(lx0+loct+17)=0.0d0
cc
cc    calculate equivalent conductances and currents for
cc    meyer"s capacitors
cc
      call intgr8(gcgs,ceqgs,cgs,loct+12)
      call intgr8(gcgd,ceqgd,cgd,loct+14)
      call intgr8(gcgb,ceqgb,cgb,loct+16)
      ceqgs=ceqgs-gcgs*vgs+ag(1)*qgs(lx0+loct)
      ceqgd=ceqgd-gcgd*vgd+ag(1)*qgd(lx0+loct)
      ceqgb=ceqgb-gcgb*vgb+ag(1)*qgb(lx0+loct)
      if (initf.ne.5) go to 870
      cqgs(lx1+loct)=cqgs(lx0+loct)
      cqgd(lx1+loct)=cqgd(lx0+loct)
      cqgb(lx1+loct)=cqgb(lx0+loct)
      go to 870
c
c.. bulk and channel charge (plus overlaps)
c
  750 if (mode.ne.1) go to 755
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 755
      if (initf.eq.4) go to 755
      go to 850
  755 if (devmod.eq.-1.0d0) go to 760
      call moscap(vgd,vgs,vgb,covlgd,covlgs,covlgb,
     1   capbd,capbs,cggb,cgdb,cgsb,cbgb,cbdb,cbsb,
     2   gcggb,gcgdb,gcgsb,gcbgb,gcbdb,gcbsb,
     3   gcdgb,gcddb,gcdsb,gcsgb,gcsdb,gcssb,
     4   qgate,qchan,qbulk,qdrn,qsrc)
      go to 780
  760 call moscap(vgs,vgd,vgb,covlgs,covlgd,covlgb,
     1   capbs,capbd,cggb,cgsb,cgdb,cbgb,cbsb,cbdb,
     2   gcggb,gcgsb,gcgdb,gcbgb,gcbsb,gcbdb,
     3   gcsgb,gcssb,gcsdb,gcdgb,gcdsb,gcddb,
     4   qgate,qchan,qbulk,qsrc,qdrn)
  780 if (ibypas.eq.1) go to 860
      qg(lx0+loct)=qgate
      qd(lx0+loct)=qdrn-qbd(lx0+loct)
      qb(lx0+loct)=qbulk+qbd(lx0+loct)+qbs(lx0+loct)
c
c  store small-signal parameters
c
  790 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 850
      if (initf.ne.4) go to 800
      value(lx0+loct+18)=cggb
      value(lx0+loct+19)=cgdb
      value(lx0+loct+20)=cgsb
      value(lx0+loct+21)=cbgb
      value(lx0+loct+22)=cbdb
      value(lx0+loct+23)=cbsb
      value(lx0+loct+24)=capbd
      value(lx0+loct+26)=capbs
      go to 1000
c
c  transient analysis
c
  800 if (initf.ne.5) go to 810
      qb(lx1+loct)=qb(lx0+loct)
      qg(lx1+loct)=qg(lx0+loct)
      qd(lx1+loct)=qd(lx0+loct)
c.. integrate qb
  810 call intgr8(geq,ceq,0.0d0,loct+12)
c.. integrate qg
      call intgr8(geq,ceq,0.0d0,loct+14)
c.. integrate qd
      call intgr8(geq,ceq,0.0d0,loct+16)
      go to 860
c
c     initialize to zero charge conductances and current
c
  850 ceqqg=0.0d0
      ceqqb=0.0d0
      ceqqd=0.0d0
      gcdgb=0.0d0
      gcddb=0.0d0
      gcdsb=0.0d0
      gcsgb=0.0d0
      gcsdb=0.0d0
      gcssb=0.0d0
      gcggb=0.0d0
      gcgdb=0.0d0
      gcgsb=0.0d0
      gcbgb=0.0d0
      gcbdb=0.0d0
      gcbsb=0.0d0
      go to 900
c
c     evaluate equivalent charge currents
c
  860 cgate=cqg(lx0+loct)
      cqbulk=cqb(lx0+loct)
      cqdrn=cqd(lx0+loct)
      ceqqg=cgate-gcggb*vgb+gcgdb*vbd+gcgsb*vbs
      ceqqb=cqbulk-gcbgb*vgb+gcbdb*vbd+gcbsb*vbs
      ceqqd=cqdrn-gcdgb*vgb+gcddb*vbd+gcdsb*vbs
      if (initf.ne.5) go to 900
      cqb(lx1+loct)=cqb(lx0+loct)
      cqg(lx1+loct)=cqg(lx0+loct)
      cqd(lx1+loct)=cqd(lx0+loct)
      go to 900
c
cc   do the mapping from meyer"s capacitor model into the charge
cc   oriented model
cc
  870 ceqqg=ceqgs+ceqgb+ceqgd
      ceqqb=-ceqgb
      ceqqd=-ceqgd
      gcbdb=0.0d0
      gcbsb=0.0d0
      gcdsb=0.0d0
      gcsdb=0.0d0
      gcgdb=-gcgd
      gcgsb=-gcgs
      gcbgb=-gcgb
      gcdgb=-gcgd
      gcsgb=-gcgs
      gcssb=gcgs
      gcddb=gcgd
      gcggb=gcgd+gcgs+gcgb
c
c     store charge storage info for meyer's cap in lx table
c
      cgbo(lx0+loct)=cgb
      cgso(lx0+loct)=cgs
      cgdo(lx0+loct)=cgd
c
c  load current vector
c
  900 ceqbs=type*(cbs-(gbs-gmin)*vbs)
      ceqbd=type*(cbd-(gbd-gmin)*vbd)
      ceqqg=type*ceqqg
      ceqqb=type*ceqqb
      ceqqd=type*ceqqd
      xnrm=1.0d0
      xrev=0.0d0
      if (devmod.lt.0.0d0) go to 910
      cdreq=type*(cdrain-gds*vds-gm*vgs-gmbs*vbs)
      go to 920
  910 xnrm=0.0d0
      xrev=1.0d0
      cdreq=-type*(cdrain-gds*(-vds)-gm*vgd-gmbs*vbd)
  920 value(lvn+node2)=value(lvn+node2)-ceqqg
      value(lvn+node4)=value(lvn+node4)-ceqbs-ceqbd-ceqqb
      value(lvn+node5)=value(lvn+node5)-cdreq+ceqbd-ceqqd
      value(lvn+node6)=value(lvn+node6)+cdreq+ceqbs
     1   +ceqqg+ceqqb+ceqqd
c
c  load y matrix
c
      locy=lvn+nodplc(loc+27)
      value(locy)=value(locy)+gdpr
      locy=lvn+nodplc(loc+28)
      value(locy)=value(locy)+gcggb
      locy=lvn+nodplc(loc+29)
      value(locy)=value(locy)+gspr
      locy=lvn+nodplc(loc+30)
      value(locy)=value(locy)+gbd+gbs-gcbgb-gcbdb-gcbsb
      locy=lvn+nodplc(loc+31)
      value(locy)=value(locy)+gdpr+gds+gbd+xrev*(gm+gmbs)+gcddb
      locy=lvn+nodplc(loc+32)
      value(locy)=value(locy)+gspr+gds+gbs+xnrm*(gm+gmbs)+gcssb
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)-gdpr
      locy=lvn+nodplc(loc+11)
      value(locy)=value(locy)-gcggb-gcgdb-gcgsb
      locy=lvn+nodplc(loc+12)
      value(locy)=value(locy)+gcgdb
      locy=lvn+nodplc(loc+13)
      value(locy)=value(locy)+gcgsb
      locy=lvn+nodplc(loc+14)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+15)
      value(locy)=value(locy)+gcbgb
      locy=lvn+nodplc(loc+16)
      value(locy)=value(locy)-gbd+gcbdb
      locy=lvn+nodplc(loc+17)
      value(locy)=value(locy)-gbs+gcbsb
      locy=lvn+nodplc(loc+18)
      value(locy)=value(locy)-gdpr
      locy=lvn+nodplc(loc+19)
      value(locy)=value(locy)+(xnrm-xrev)*gm+gcdgb
      locy=lvn+nodplc(loc+20)
      value(locy)=value(locy)-gbd+(xnrm-xrev)*gmbs-
     1   gcdgb-gcddb-gcdsb
      locy=lvn+nodplc(loc+21)
      value(locy)=value(locy)-gds-xnrm*(gm+gmbs)+gcdsb
      locy=lvn+nodplc(loc+22)
      value(locy)=value(locy)-(xnrm-xrev)*gm+gcsgb
      locy=lvn+nodplc(loc+23)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+24)
      value(locy)=value(locy)-gbs-(xnrm-xrev)*gmbs-
     1   gcsgb-gcsdb-gcssb
      locy=lvn+nodplc(loc+25)
      value(locy)=value(locy)-gds-xrev*(gm+gmbs)+gcsdb
 1000 loc=nodplc(loc)
      go to 10
      end
