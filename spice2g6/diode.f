      subroutine diode
      implicit double precision (a-h,o-z)
c
c     this routine processes diodes for dc and transient analyses.
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
      dimension vdo(1),cdo(1),gdo(1),qd(1),cqd(1)
      equivalence (vdo(1),value(1)),(cdo(1),value(2)),
     1   (gdo(1),value(3)),(qd(1),value(4)),(cqd(1),value(5))
c
c
      loc=locate(11)
   10 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) return
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+5)
      ioff=nodplc(loc+6)
      locm=nodplc(locm+1)
      loct=nodplc(loc+11)
c
c  dc model parameters
c
      area=value(locv+1)
      csat=value(locm+1)*area
      gspr=value(locm+2)*area
      vte=value(locm+3)*vt
      bv=value(locm+13)
      vcrit=value(locm+18)
c
c  initialization
c
      icheck=1
      go to (100,20,30,50,60,70),initf
   20 if(mode.ne.1.or.modedc.ne.2.or.nosolv.eq.0) go to 25
      vd=value(locv+2)
      go to 300
   25 if(ioff.ne.0) go to 40
      vd=vcrit
      go to 300
   30 if (ioff.eq.0) go to 100
   40 vd=0.0d0
      go to 300
   50 vd=vdo(lx0+loct)
      go to 300
   60 vd=vdo(lx1+loct)
      go to 300
   70 xfact=delta/delold(2)
      vdo(lx0+loct)=vdo(lx1+loct)
      vd=(1.0d0+xfact)*vdo(lx1+loct)-xfact*vdo(lx2+loct)
      cdo(lx0+loct)=cdo(lx1+loct)
      gdo(lx0+loct)=gdo(lx1+loct)
      go to 110
c
c  compute new nonlinear branch voltage
c
  100 vd=value(lvnim1+node3)-value(lvnim1+node2)
  110 delvd=vd-vdo(lx0+loct)
      cdhat=cdo(lx0+loct)+gdo(lx0+loct)*delvd
c
c  bypass if solution has not changed
c
      if (initf.eq.6) go to 200
      tol=reltol*dmax1(dabs(vd),dabs(vdo(lx0+loct)))+vntol
      if (dabs(delvd).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cdhat),dabs(cdo(lx0+loct)))+abstol
      if (dabs(cdhat-cdo(lx0+loct)).ge.tol) go to 200
      vd=vdo(lx0+loct)
      cd=cdo(lx0+loct)
      gd=gdo(lx0+loct)
      go to 800
c
c  limit new junction voltage
c
  200 vlim=vte+vte
      if(bv.eq.0.0d0) go to 205
      if (vd.lt.dmin1(0.0d0,-bv+10.0d0*vte)) go to 210
  205 call pnjlim(vd,vdo(lx0+loct),vte,vcrit,icheck)
      go to 300
  210 vdtemp=-(vd+bv)
      call pnjlim(vdtemp,-(vdo(lx0+loct)+bv),vte,vcrit,icheck)
      vd=-(vdtemp+bv)
c
c  compute dc current and derivitives
c
  300 if (vd.lt.-5.0d0*vte) go to 310
      evd=dexp(vd/vte)
      cd=csat*(evd-1.0d0)+gmin*vd
      gd=csat*evd/vte+gmin
      go to 330
  310 if(bv.eq.0.0d0) go to 315
      if(vd.lt.-bv) go to 320
  315 gd=-csat/vd+gmin
      cd=gd*vd
      go to 330
  320 evrev=dexp(-(bv+vd)/vt)
      cd=-csat*(evrev-1.0d0+bv/vt)
      gd=csat*evrev/vt
  330 if (mode.ne.1) go to 500
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 500
      if (initf.eq.4) go to 500
      go to 700
c
c  charge storage elements
c
  500 tau=value(locm+4)
      czero=value(locm+5)*area
      pb=value(locm+6)
      xm=value(locm+7)
      fcpb=value(locm+12)
      if (vd.ge.fcpb) go to 510
      arg=1.0d0-vd/pb
      sarg=dexp(-xm*dlog(arg))
      qd(lx0+loct)=tau*cd+pb*czero*(1.0d0-arg*sarg)/(1.0d0-xm)
      capd=tau*gd+czero*sarg
      go to 520
  510 f1=value(locm+15)
      f2=value(locm+16)
      f3=value(locm+17)
      czof2=czero/f2
      qd(lx0+loct)=tau*cd+czero*f1+czof2*(f3*(vd-fcpb)
     1   +(xm/(pb+pb))*(vd*vd-fcpb*fcpb))
      capd=tau*gd+czof2*(f3+xm*vd/pb)
c
c  store small-signal parameters
c
  520 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 700
      if (initf.ne.4) go to 600
      value(lx0+loct+4)=capd
      go to 1000
c
c  transient analysis
c
  600 if (initf.ne.5) go to 610
      qd(lx1+loct)=qd(lx0+loct)
  610 call intgr8(geq,ceq,capd,loct+3)
      gd=gd+geq
      cd=cd+cqd(lx0+loct)
      if (initf.ne.5) go to 700
      cqd(lx1+loct)=cqd(lx0+loct)
c
c  check convergence
c
  700 if (initf.ne.3) go to 710
      if (ioff.eq.0) go to 710
      go to 750
  710 if (icheck.eq.1) go to 720
      tol=reltol*dmax1(dabs(cdhat),dabs(cd))+abstol
      if (dabs(cdhat-cd).le.tol) go to 750
  720 noncon=noncon+1
  750 vdo(lx0+loct)=vd
      cdo(lx0+loct)=cd
      gdo(lx0+loct)=gd
c
c  load current vector
c
  800 cdeq=cd-gd*vd
      value(lvn+node2)=value(lvn+node2)+cdeq
      value(lvn+node3)=value(lvn+node3)-cdeq
c
c  load matrix
c
      locy=lvn+nodplc(loc+13)
      value(locy)=value(locy)+gspr
      locy=lvn+nodplc(loc+14)
      value(locy)=value(locy)+gd
      locy=lvn+nodplc(loc+15)
      value(locy)=value(locy)+gd+gspr
      locy=lvn+nodplc(loc+7)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+8)
      value(locy)=value(locy)-gd
      locy=lvn+nodplc(loc+9)
      value(locy)=value(locy)-gspr
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)-gd
 1000 loc=nodplc(loc)
      go to 10
      end
