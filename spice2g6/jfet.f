      subroutine jfet
      implicit double precision (a-h,o-z)
c
c     this routine processes jfets for dc and transient analyses.
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension vgso(1),vgdo(1),cgo(1),cdo(1),cgdo(1),gmo(1),gdso(1),
     1   ggso(1),ggdo(1),qgs(1),cqgs(1),qgd(1),cqgd(1)
      equivalence (vgso(1),value( 1)),(vgdo(1),value( 2)),
     1            (cgo (1),value( 3)),(cdo (1),value( 4)),
     2            (cgdo(1),value( 5)),(gmo (1),value( 6)),
     3            (gdso(1),value( 7)),(ggso(1),value( 8)),
     4            (ggdo(1),value( 9)),(qgs (1),value(10)),
     5            (cqgs(1),value(11)),(qgd (1),value(12)),
     6            (cqgd(1),value(13))
c
c
      loc=locate(13)
   10 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) return
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      locm=nodplc(loc+7)
      ioff=nodplc(loc+8)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=nodplc(loc+19)
c
c  dc model parameters
c
      area=value(locv+1)
      vto=value(locm+1)
      beta=value(locm+2)*area
      xlamb=value(locm+3)
      gdpr=value(locm+4)*area
      gspr=value(locm+5)*area
      csat=value(locm+9)*area
      vcrit=value(locm+16)
c
c  initialization
c
      icheck=1
      go to (100,20,30,50,60,70), initf
   20 if(mode.ne.1.or.modedc.ne.2.or.nosolv.eq.0) go to 25
      vds=type*value(locv+2)
      vgs=type*value(locv+3)
      vgd=vgs-vds
      go to 300
   25 if(ioff.ne.0) go to 40
      vgs=-1.0d0
      vgd=-1.0d0
      go to 300
   30 if (ioff.eq.0) go to 100
   40 vgs=0.0d0
      vgd=0.0d0
      go to 300
   50 vgs=vgso(lx0+loct)
      vgd=vgdo(lx0+loct)
      go to 300
   60 vgs=vgso(lx1+loct)
      vgd=vgdo(lx1+loct)
      go to 300
   70 xfact=delta/delold(2)
      vgso(lx0+loct)=vgso(lx1+loct)
      vgs=(1.0d0+xfact)*vgso(lx1+loct)-xfact*vgso(lx2+loct)
      vgdo(lx0+loct)=vgdo(lx1+loct)
      vgd=(1.0d0+xfact)*vgdo(lx1+loct)-xfact*vgdo(lx2+loct)
      cgo(lx0+loct)=cgo(lx1+loct)
      cdo(lx0+loct)=cdo(lx1+loct)
      cgdo(lx0+loct)=cgdo(lx1+loct)
      gmo(lx0+loct)=gmo(lx1+loct)
      gdso(lx0+loct)=gdso(lx1+loct)
      ggso(lx0+loct)=ggso(lx1+loct)
      ggdo(lx0+loct)=ggdo(lx1+loct)
      go to 110
c
c  compute new nonlinear branch voltages
c
  100 vgs=type*(value(lvnim1+node2)-value(lvnim1+node5))
      vgd=type*(value(lvnim1+node2)-value(lvnim1+node4))
  110 delvgs=vgs-vgso(lx0+loct)
      delvgd=vgd-vgdo(lx0+loct)
      delvds=delvgs-delvgd
      cghat=cgo(lx0+loct)+ggdo(lx0+loct)*delvgd+ggso(lx0+loct)*delvgs
      cdhat=cdo(lx0+loct)+gmo(lx0+loct)*delvgs+gdso(lx0+loct)*delvds
     1   -ggdo(lx0+loct)*delvgd
c
c  bypass if solution has not changed
c
      if (initf.eq.6) go to 200
      tol=reltol*dmax1(dabs(vgs),dabs(vgso(lx0+loct)))+vntol
      if (dabs(delvgs).ge.tol) go to 200
      tol=reltol*dmax1(dabs(vgd),dabs(vgdo(lx0+loct)))+vntol
      if (dabs(delvgd).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cghat),dabs(cgo(lx0+loct)))+abstol
      if (dabs(cghat-cgo(lx0+loct)).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cdhat),dabs(cdo(lx0+loct)))+abstol
      if (dabs(cdhat-cdo(lx0+loct)).ge.tol) go to 200
      vgs=vgso(lx0+loct)
      vgd=vgdo(lx0+loct)
      vds=vgs-vgd
      cg=cgo(lx0+loct)
      cd=cdo(lx0+loct)
      cgd=cgdo(lx0+loct)
      gm=gmo(lx0+loct)
      gds=gdso(lx0+loct)
      ggs=ggso(lx0+loct)
      ggd=ggdo(lx0+loct)
      go to 900
c
c  limit nonlinear branch voltages
c
  200 ichk1=1
      call pnjlim(vgs,vgso(lx0+loct),vt,vcrit,icheck)
      call pnjlim(vgd,vgdo(lx0+loct),vt,vcrit,ichk1)
      if (ichk1.eq.1) icheck=1
      call fetlim(vgs,vgso(lx0+loct),vto)
      call fetlim(vgd,vgdo(lx0+loct),vto)
c
c  determine dc current and derivatives
c
  300 vds=vgs-vgd
      if (vgs.gt.-5.0d0*vt) go to 310
      ggs=-csat/vgs+gmin
      cg=ggs*vgs
      go to 320
  310 evgs=dexp(vgs/vt)
      ggs=csat*evgs/vt+gmin
      cg=csat*(evgs-1.0d0)+gmin*vgs
  320 if (vgd.gt.-5.0d0*vt) go to 330
      ggd=-csat/vgd+gmin
      cgd=ggd*vgd
      go to 340
  330 evgd=dexp(vgd/vt)
      ggd=csat*evgd/vt+gmin
      cgd=csat*(evgd-1.0d0)+gmin*vgd
  340 cg=cg+cgd
c
c  compute drain current and derivitives for normal mode
c
  400 if (vds.lt.0.0d0) go to 450
      vgst=vgs-vto
c
c  normal mode, cutoff region
c
      if (vgst.gt.0.0d0) go to 410
      cdrain=0.0d0
      gm=0.0d0
      gds=0.0d0
      go to 490
c
c  normal mode, saturation region
c
  410 betap=beta*(1.0d0+xlamb*vds)
      twob=betap+betap
      if (vgst.gt.vds) go to 420
      cdrain=betap*vgst*vgst
      gm=twob*vgst
      gds=xlamb*beta*vgst*vgst
      go to 490
c
c  normal mode, linear region
c
  420 cdrain=betap*vds*(vgst+vgst-vds)
      gm=twob*vds
      gds=twob*(vgst-vds)+xlamb*beta*vds*(vgst+vgst-vds)
      go to 490
c
c  compute drain current and derivitives for inverse mode
c
  450 vgdt=vgd-vto
c
c  inverse mode, cutoff region
c
      if (vgdt.gt.0.0d0) go to 460
      cdrain=0.0d0
      gm=0.0d0
      gds=0.0d0
      go to 490
c
c  inverse mode, saturation region
c
  460 betap=beta*(1.0d0-xlamb*vds)
      twob=betap+betap
      if (vgdt.gt.-vds) go to 470
      cdrain=-betap*vgdt*vgdt
      gm=-twob*vgdt
      gds=xlamb*beta*vgdt*vgdt-gm
      go to 490
c
c  inverse mode, linear region
c
  470 cdrain=betap*vds*(vgdt+vgdt+vds)
      gm=twob*vds
      gds=twob*vgdt-xlamb*beta*vds*(vgdt+vgdt+vds)
c
c  compute equivalent drain current source
c
  490 cd=cdrain-cgd
      if (mode.ne.1) go to 500
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 500
      if (initf.eq.4) go to 500
      go to 700
c
c  charge storage elements
c
  500 czgs=value(locm+6)*area
      czgd=value(locm+7)*area
      phib=value(locm+8)
      twop=phib+phib
      fcpb=value(locm+12)
      fcpb2=fcpb*fcpb
      f1=value(locm+13)
      f2=value(locm+14)
      f3=value(locm+15)
      czgsf2=czgs/f2
      czgdf2=czgd/f2
      if (vgs.ge.fcpb) go to 510
      sarg=dsqrt(1.0d0-vgs/phib)
      qgs(lx0+loct)=twop*czgs*(1.0d0-sarg)
      capgs=czgs/sarg
      go to 520
  510 qgs(lx0+loct)=czgs*f1+czgsf2*(f3*(vgs-fcpb)
     1   +(vgs*vgs-fcpb2)/(twop+twop))
      capgs=czgsf2*(f3+vgs/twop)
  520 if (vgd.ge.fcpb) go to 530
      sarg=dsqrt(1.0d0-vgd/phib)
      qgd(lx0+loct)=twop*czgd*(1.0d0-sarg)
      capgd=czgd/sarg
      go to 560
  530 qgd(lx0+loct)=czgd*f1+czgdf2*(f3*(vgd-fcpb)
     1   +(vgd*vgd-fcpb2)/(twop+twop))
      capgd=czgdf2*(f3+vgd/twop)
c
c  store small-signal parameters
c
  560 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 700
      if (initf.ne.4) go to 600
      value(lx0+loct+9)=capgs
      value(lx0+loct+11)=capgd
      go to 1000
c
c  transient analysis
c
  600 if (initf.ne.5) go to 610
      qgs(lx1+loct)=qgs(lx0+loct)
      qgd(lx1+loct)=qgd(lx0+loct)
  610 call intgr8(geq,ceq,capgs,loct+9)
      ggs=ggs+geq
      cg=cg+cqgs(lx0+loct)
      call intgr8(geq,ceq,capgd,loct+11)
      ggd=ggd+geq
      cg=cg+cqgd(lx0+loct)
      cd=cd-cqgd(lx0+loct)
      cgd=cgd+cqgd(lx0+loct)
      if (initf.ne.5) go to 700
      cqgs(lx1+loct)=cqgs(lx0+loct)
      cqgd(lx1+loct)=cqgd(lx0+loct)
c
c  check convergence
c
  700 if (initf.ne.3) go to 710
      if (ioff.eq.0) go to 710
      go to 750
  710 if (icheck.eq.1) go to 720
      tol=reltol*dmax1(dabs(cghat),dabs(cg))+abstol
      if (dabs(cghat-cg).ge.tol) go to 720
      tol=reltol*dmax1(dabs(cdhat),dabs(cd))+abstol
      if (dabs(cdhat-cd).le.tol) go to 750
  720 noncon=noncon+1
  750 vgso(lx0+loct)=vgs
      vgdo(lx0+loct)=vgd
      cgo(lx0+loct)=cg
      cdo(lx0+loct)=cd
      cgdo(lx0+loct)=cgd
      gmo(lx0+loct)=gm
      gdso(lx0+loct)=gds
      ggso(lx0+loct)=ggs
      ggdo(lx0+loct)=ggd
c
c  load current vector
c
  900 ceqgd=type*(cgd-ggd*vgd)
      ceqgs=type*((cg-cgd)-ggs*vgs)
      cdreq=type*((cd+cgd)-gds*vds-gm*vgs)
      value(lvn+node2)=value(lvn+node2)-ceqgs-ceqgd
      value(lvn+node4)=value(lvn+node4)-cdreq+ceqgd
      value(lvn+node5)=value(lvn+node5)+cdreq+ceqgs
c
c  load y matrix
c
      locy=lvn+nodplc(loc+20)
      value(locy)=value(locy)+gdpr
      locy=lvn+nodplc(loc+21)
      value(locy)=value(locy)+ggd+ggs
      locy=lvn+nodplc(loc+22)
      value(locy)=value(locy)+gspr
      locy=lvn+nodplc(loc+23)
      value(locy)=value(locy)+gdpr+gds+ggd
      locy=lvn+nodplc(loc+24)
      value(locy)=value(locy)+gspr+gds+gm+ggs
      locy=lvn+nodplc(loc+9)
      value(locy)=value(locy)-gdpr
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)-ggd
      locy=lvn+nodplc(loc+11)
      value(locy)=value(locy)-ggs
      locy=lvn+nodplc(loc+12)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+13)
      value(locy)=value(locy)-gdpr
      locy=lvn+nodplc(loc+14)
      value(locy)=value(locy)+gm-ggd
      locy=lvn+nodplc(loc+15)
      value(locy)=value(locy)-gds-gm
      locy=lvn+nodplc(loc+16)
      value(locy)=value(locy)-ggs-gm
      locy=lvn+nodplc(loc+17)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+18)
      value(locy)=value(locy)-gds
 1000 loc=nodplc(loc)
      go to 10
      end
