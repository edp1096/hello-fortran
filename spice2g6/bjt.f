      subroutine bjt
      implicit double precision (a-h,o-z)
c
c     this routine processes bjts for dc and transient analyses.
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
      dimension vbeo(1),vbco(1),cco(1),cbo(1),gpio(1),gmuo(1),gmo(1),
     1   goo(1),qbe(1),cqbe(1),qbc(1),cqbc(1),qcs(1),cqcs(1),qbx(1),
     2   cqbx(1),gxo(1),cexbc(1),geqcbo(1)
      equivalence (vbeo(1),value(1)),(vbco(1),value(2)),
     1   (cco(1),value(3)),(cbo(1),value(4)),(gpio(1),value(5)),
     2   (gmuo(1),value(6)),(gmo(1),value(7)),(goo(1),value(8)),
     3   (qbe(1),value(9)),(cqbe(1),value(10)),(qbc(1),value(11)),
     4   (cqbc(1),value(12)),(qcs(1),value(13)),(cqcs(1),value(14)),
     5   (qbx(1),value(15)),(cqbx(1),value(16)),(gxo(1),value(17)),
     6   (cexbc(1),value(18)),(geqcbo(1),value(19))
c
c
      loc=locate(12)
   10 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) return
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      node7=nodplc(loc+30)
      locm=nodplc(loc+8)
      ioff=nodplc(loc+9)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=nodplc(loc+22)
      gccs=0.0d0
      ceqcs=0.0d0
      geqbx=0.0d0
      ceqbx=0.0d0
      geqcb=0.0d0
c
c  dc model paramters
c
      area=value(locv+1)
      bfm=value(locm+2)
      brm=value(locm+8)
      csat=value(locm+1)*area
      rbpr=value(locm+18)/area
      rbpi=value(locm+16)/area-rbpr
      gcpr=value(locm+20)*area
      gepr=value(locm+19)*area
      ova=value(locm+4)
      ovb=value(locm+10)
      oik=value(locm+5)/area
      c2=value(locm+6)*area
      vte=value(locm+7)*vt
      oikr=value(locm+11)/area
      c4=value(locm+12)*area
      vtc=value(locm+13)*vt
      vcrit=value(locm+54)
      td=value(locm+28)
      xjrb=value(locm+17)*area
c
c  initialization
c
      icheck=1
      go to (100,20,30,50,60,70),initf
   20 if(mode.ne.1.or.modedc.ne.2.or.nosolv.eq.0) go to 25
      vbe=type*value(locv+2)
      vce=type*value(locv+3)
      vbc=vbe-vce
      vbx=vbc
      vcs=0.0d0
      go to 300
   25 if(ioff.ne.0) go to 40
      vbe=vcrit
      vbc=0.0d0
      go to 300
   30 if (ioff.eq.0) go to 100
   40 vbe=0.0d0
      vbc=0.0d0
      go to 300
   50 vbe=vbeo(lx0+loct)
      vbc=vbco(lx0+loct)
      vbx=type*(value(lvnim1+node2)-value(lvnim1+node4))
      vcs=type*(value(lvnim1+node7)-value(lvnim1+node4))
      go to 300
   60 vbe=vbeo(lx1+loct)
      vbc=vbco(lx1+loct)
      vbx=type*(value(lvnim1+node2)-value(lvnim1+node4))
      vcs=type*(value(lvnim1+node7)-value(lvnim1+node4))
      if(mode.ne.2.or.nosolv.eq.0) go to 300
      vbx=type*(value(locv+2)-value(locv+3))
      vcs=0.0d0
      go to 300
   70 xfact=delta/delold(2)
      vbeo(lx0+loct)=vbeo(lx1+loct)
      vbe=(1.0d0+xfact)*vbeo(lx1+loct)-xfact*vbeo(lx2+loct)
      vbco(lx0+loct)=vbco(lx1+loct)
      vbc=(1.0d0+xfact)*vbco(lx1+loct)-xfact*vbco(lx2+loct)
      cco(lx0+loct)=cco(lx1+loct)
      cbo(lx0+loct)=cbo(lx1+loct)
      gpio(lx0+loct)=gpio(lx1+loct)
      gmuo(lx0+loct)=gmuo(lx1+loct)
      gmo(lx0+loct)=gmo(lx1+loct)
      goo(lx0+loct)=goo(lx1+loct)
      gxo(lx0+loct)=gxo(lx1+loct)
      go to 110
c
c  compute new nonlinear branch voltages
c
  100 vbe=type*(value(lvnim1+node5)-value(lvnim1+node6))
      vbc=type*(value(lvnim1+node5)-value(lvnim1+node4))
  110 delvbe=vbe-vbeo(lx0+loct)
      delvbc=vbc-vbco(lx0+loct)
      vbx=type*(value(lvnim1+node2)-value(lvnim1+node4))
      vcs=type*(value(lvnim1+node7)-value(lvnim1+node4))
      cchat=cco(lx0+loct)+(gmo(lx0+loct)+goo(lx0+loct))*delvbe
     1   -(goo(lx0+loct)+gmuo(lx0+loct))*delvbc
      cbhat=cbo(lx0+loct)+gpio(lx0+loct)*delvbe+gmuo(lx0+loct)*delvbc
c
c   bypass if solution has not changed
c
      if (initf.eq.6) go to 200
      tol=reltol*dmax1(dabs(vbe),dabs(vbeo(lx0+loct)))+vntol
      if (dabs(delvbe).ge.tol) go to 200
      tol=reltol*dmax1(dabs(vbc),dabs(vbco(lx0+loct)))+vntol
      if (dabs(delvbc).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cchat),dabs(cco(lx0+loct)))+abstol
      if (dabs(cchat-cco(lx0+loct)).ge.tol) go to 200
      tol=reltol*dmax1(dabs(cbhat),dabs(cbo(lx0+loct)))+abstol
      if (dabs(cbhat-cbo(lx0+loct)).ge.tol) go to 200
      vbe=vbeo(lx0+loct)
      vbc=vbco(lx0+loct)
      cc=cco(lx0+loct)
      cb=cbo(lx0+loct)
      gpi=gpio(lx0+loct)
      gmu=gmuo(lx0+loct)
      gm=gmo(lx0+loct)
      go=goo(lx0+loct)
      gx=gxo(lx0+loct)
      geqcb=geqcbo(lx0+loct)
      if (mode.ne.1) go to 800
      go to 900
c
c  limit nonlinear branch voltages
c
  200 ichk1=1
      call pnjlim(vbe,vbeo(lx0+loct),vt,vcrit,icheck)
      call pnjlim(vbc,vbco(lx0+loct),vt,vcrit,ichk1)
      if (ichk1.eq.1) icheck=1
c
c  determine dc current and derivitives
c
  300 vtn=vt*value(locm+3)
      if(vbe.le.-5.0d0*vtn) go to 320
      evbe=dexp(vbe/vtn)
      cbe=csat*(evbe-1.0d0)+gmin*vbe
      gbe=csat*evbe/vtn+gmin
      if (c2.ne.0.0d0) go to 310
      cben=0.0d0
      gben=0.0d0
      go to 350
  310 evben=dexp(vbe/vte)
      cben=c2*(evben-1.0d0)
      gben=c2*evben/vte
      go to 350
  320 gbe=-csat/vbe+gmin
      cbe=gbe*vbe
      gben=-c2/vbe
      cben=gben*vbe
  350 vtn=vt*value(locm+9)
      if(vbc.le.-5.0d0*vtn) go to 370
      evbc=dexp(vbc/vtn)
      cbc=csat*(evbc-1.0d0)+gmin*vbc
      gbc=csat*evbc/vtn+gmin
      if (c4.ne.0.0d0) go to 360
      cbcn=0.0d0
      gbcn=0.0d0
      go to 400
  360 evbcn=dexp(vbc/vtc)
      cbcn=c4*(evbcn-1.0d0)
      gbcn=c4*evbcn/vtc
      go to 400
  370 gbc=-csat/vbc+gmin
      cbc=gbc*vbc
      gbcn=-c4/vbc
      cbcn=gbcn*vbc
c
c  determine base charge terms
c
  400 q1=1.0d0/(1.0d0-ova*vbc-ovb*vbe)
      if (oik.ne.0.0d0) go to 405
      if (oikr.ne.0.0d0) go to 405
      qb=q1
      dqbdve=q1*qb*ovb
      dqbdvc=q1*qb*ova
      go to 410
  405 q2=oik*cbe+oikr*cbc
      arg=dmax1(0.0d0,1.0d0+4.0d0*q2)
      sqarg=1.0d0
      if(arg.ne.0.0d0) sqarg=dsqrt(arg)
      qb=q1*(1.0d0+sqarg)/2.0d0
      dqbdve=q1*(qb*ovb+oik*gbe/sqarg)
      dqbdvc=q1*(qb*ova+oikr*gbc/sqarg)
c
c  weil's approx. for excess phase applied with backward-
c  euler integration
c
  410 cc=0.0d0
      cex=cbe
      gex=gbe
      if(mode.eq.1) go to 420
      if(td.eq.0.0d0) go to 420
      arg1=delta/td
      arg2=3.0d0*arg1
      arg1=arg2*arg1
      denom=1.0d0+arg1+arg2
      arg3=arg1/denom
      if(initf.ne.5) go to 411
      cexbc(lx1+loct)=cbe/qb
      cexbc(lx2+loct)=cexbc(lx1+loct)
  411 cc=(cexbc(lx1+loct)*(1.0d0+delta/delold(2)+arg2)
     1  -cexbc(lx2+loct)*delta/delold(2))/denom
      cex=cbe*arg3
      gex=gbe*arg3
      cexbc(lx0+loct)=cc+cex/qb
c
c  determine dc incremental conductances
c
  420 cc=cc+(cex-cbc)/qb-cbc/brm-cbcn
      cb=cbe/bfm+cben+cbc/brm+cbcn
      gx=rbpr+rbpi/qb
      if(xjrb.eq.0.0d0) go to 430
      arg1=dmax1(cb/xjrb,1.0d-9)
      arg2=(-1.0d0+dsqrt(1.0d0+14.59025d0*arg1))/2.4317d0/dsqrt(arg1)
      arg1=dtan(arg2)
      gx=rbpr+3.0d0*rbpi*(arg1-arg2)/arg2/arg1/arg1
  430 if(gx.ne.0.0d0) gx=1.0d0/gx
      gpi=gbe/bfm+gben
      gmu=gbc/brm+gbcn
      go=(gbc+(cex-cbc)*dqbdvc/qb)/qb
      gm=(gex-(cex-cbc)*dqbdve/qb)/qb-go
      if (mode.ne.1) go to 500
      if ((modedc.eq.2).and.(nosolv.ne.0)) go to 500
      if (initf.eq.4) go to 500
      go to 700
c
c  charge storage elements
c
  500 tf=value(locm+24)
      tr=value(locm+33)
      czbe=value(locm+21)*area
      pe=value(locm+22)
      xme=value(locm+23)
      cdis=value(locm+32)
      ctot=value(locm+29)*area
      czbc=ctot*cdis
      czbx=ctot-czbc
      pc=value(locm+30)
      xmc=value(locm+31)
      fcpe=value(locm+46)
      czcs=value(locm+38)*area
      ps=value(locm+39)
      xms=value(locm+40)
      xtf=value(locm+25)
      ovtf=value(locm+26)
      xjtf=value(locm+27)*area
      if(tf.eq.0.0d0) go to 505
      if(vbe.le.0.0d0) go to 505
      argtf=0.0d0
      arg2=0.0d0
      arg3=0.0d0
      if(xtf.eq.0.0d0) go to 504
      argtf=xtf
      if(ovtf.ne.0.0d0) argtf=argtf*dexp(vbc*ovtf)
      arg2=argtf
      if(xjtf.eq.0.0d0) go to 503
      temp=cbe/(cbe+xjtf)
      argtf=argtf*temp*temp
      arg2=argtf*(3.0d0-temp-temp)
  503 arg3=cbe*argtf*ovtf
  504 cbe=cbe*(1.0d0+argtf)/qb
      gbe=(gbe*(1.0d0+arg2)-cbe*dqbdve)/qb
      geqcb=tf*(arg3-cbe*dqbdvc)/qb
  505 if (vbe.ge.fcpe) go to 510
      arg=1.0d0-vbe/pe
      sarg=dexp(-xme*dlog(arg))
      qbe(lx0+loct)=tf*cbe+pe*czbe*(1.0d0-arg*sarg)/(1.0d0-xme)
      capbe=tf*gbe+czbe*sarg
      go to 520
  510 f1=value(locm+47)
      f2=value(locm+48)
      f3=value(locm+49)
      czbef2=czbe/f2
      qbe(lx0+loct)=tf*cbe+czbe*f1+czbef2*(f3*(vbe-fcpe)
     1   +(xme/(pe+pe))*(vbe*vbe-fcpe*fcpe))
      capbe=tf*gbe+czbef2*(f3+xme*vbe/pe)
  520 fcpc=value(locm+50)
      f1=value(locm+51)
      f2=value(locm+52)
      f3=value(locm+53)
      if (vbc.ge.fcpc) go to 530
      arg=1.0d0-vbc/pc
      sarg=dexp(-xmc*dlog(arg))
      qbc(lx0+loct)=tr*cbc+pc*czbc*(1.0d0-arg*sarg)/(1.0d0-xmc)
      capbc=tr*gbc+czbc*sarg
      go to 540
  530 czbcf2=czbc/f2
      qbc(lx0+loct)=tr*cbc+czbc*f1+czbcf2*(f3*(vbc-fcpc)
     1   +(xmc/(pc+pc))*(vbc*vbc-fcpc*fcpc))
      capbc=tr*gbc+czbcf2*(f3+xmc*vbc/pc)
  540 if(vbx.ge.fcpc) go to 550
      arg=1.0d0-vbx/pc
      sarg=dexp(-xmc*dlog(arg))
      qbx(lx0+loct)=pc*czbx*(1.0d0-arg*sarg)/(1.0d0-xmc)
      capbx=czbx*sarg
      go to 560
  550 czbxf2=czbx/f2
      qbx(lx0+loct)=czbx*f1+czbxf2*(f3*(vbx-fcpc)+(xmc/(pc+pc))*
     1   (vbx*vbx-fcpc*fcpc))
      capbx=czbxf2*(f3+xmc*vbx/pc)
  560 if(vcs.ge.0.0d0) go to 570
      arg=1.0d0-vcs/ps
      sarg=dexp(-xms*dlog(arg))
      qcs(lx0+loct)=ps*czcs*(1.0d0-arg*sarg)/(1.0d0-xms)
      capcs=czcs*sarg
      go to 580
  570 qcs(lx0+loct)=vcs*czcs*(1.0d0+xms*vcs/(2.0d0*ps))
      capcs=czcs*(1.0d0+xms*vcs/ps)
c
c  store small-signal parameters
c
  580 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 700
      if (initf.ne.4) go to 600
      value(lx0+loct+9)=capbe
      value(lx0+loct+11)=capbc
      value(lx0+loct+13)=capcs
      value(lx0+loct+15)=capbx
      value(lx0+loct+17)=geqcb
      go to 1000
c
c  transient analysis
c
  600 if (initf.ne.5) go to 610
      qbe(lx1+loct)=qbe(lx0+loct)
      qbc(lx1+loct)=qbc(lx0+loct)
      qbx(lx1+loct)=qbx(lx0+loct)
      qcs(lx1+loct)=qcs(lx0+loct)
  610 call intgr8(geq,ceq,capbe,loct+8)
      geqcb=geqcb*ag(1)
      gpi=gpi+geq
      cb=cb+cqbe(lx0+loct)
      call intgr8(geq,ceq,capbc,loct+10)
      gmu=gmu+geq
      cb=cb+cqbc(lx0+loct)
      cc=cc-cqbc(lx0+loct)
      if (initf.ne.5) go to 700
      cqbe(lx1+loct)=cqbe(lx0+loct)
      cqbc(lx1+loct)=cqbc(lx0+loct)
c
c  check convergence
c
  700 if (initf.ne.3) go to 710
      if (ioff.eq.0) go to 710
      go to 750
  710 if (icheck.eq.1) go to 720
      tol=reltol*dmax1(dabs(cchat),dabs(cc))+abstol
      if (dabs(cchat-cc).gt.tol) go to 720
      tol=reltol*dmax1(dabs(cbhat),dabs(cb))+abstol
      if (dabs(cbhat-cb).le.tol) go to 750
  720 noncon=noncon+1
  750 vbeo(lx0+loct)=vbe
      vbco(lx0+loct)=vbc
      cco(lx0+loct)=cc
      cbo(lx0+loct)=cb
      gpio(lx0+loct)=gpi
      gmuo(lx0+loct)=gmu
      gmo(lx0+loct)=gm
      goo(lx0+loct)=go
      gxo(lx0+loct)=gx
      geqcbo(lx0+loct)=geqcb
      if (mode.eq.1) go to 900
c
c     charge storage for c-s and b-x junctions
c
  800 call intgr8(gccs,ceq,capcs,loct+12)
      ceqcs=type*(cqcs(lx0+loct)-vcs*gccs)
      call intgr8(geqbx,ceq,capbx,loct+14)
      ceqbx=type*(cqbx(lx0+loct)-vbx*geqbx)
      if (initf.ne.5) go to 900
      cqbx(lx1+loct)=cqbx(lx0+loct)
      cqcs(lx1+loct)=cqcs(lx0+loct)
c
c  load current excitation vector
c
  900 ceqbe=type*(cc+cb-vbe*(gm+go+gpi)+vbc*(go-geqcb))
      ceqbc=type*(-cc+vbe*(gm+go)-vbc*(gmu+go))
      value(lvn+node2)=value(lvn+node2)-ceqbx
      value(lvn+node4)=value(lvn+node4)+ceqcs+ceqbx+ceqbc
      value(lvn+node5)=value(lvn+node5)-ceqbe-ceqbc
      value(lvn+node6)=value(lvn+node6)+ceqbe
      value(lvn+node7)=value(lvn+node7)-ceqcs
c
c  load y matrix
c
      locy=lvn+nodplc(loc+24)
      value(locy)=value(locy)+gcpr
      locy=lvn+nodplc(loc+25)
      value(locy)=value(locy)+gx+geqbx
      locy=lvn+nodplc(loc+26)
      value(locy)=value(locy)+gepr
      locy=lvn+nodplc(loc+27)
      value(locy)=value(locy)+gmu+go+gcpr+gccs+geqbx
      locy=lvn+nodplc(loc+28)
      value(locy)=value(locy)+gx  +gpi+gmu+geqcb
      locy=lvn+nodplc(loc+29)
      value(locy)=value(locy)+gpi+gepr+gm+go
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)-gcpr
      locy=lvn+nodplc(loc+11)
      value(locy)=value(locy)-gx
      locy=lvn+nodplc(loc+12)
      value(locy)=value(locy)-gepr
      locy=lvn+nodplc(loc+13)
      value(locy)=value(locy)-gcpr
      locy=lvn+nodplc(loc+14)
      value(locy)=value(locy)-gmu+gm
      locy=lvn+nodplc(loc+15)
      value(locy)=value(locy)-gm-go
      locy=lvn+nodplc(loc+16)
      value(locy)=value(locy)-gx
      locy=lvn+nodplc(loc+17)
      value(locy)=value(locy)-gmu-geqcb
      locy=lvn+nodplc(loc+18)
      value(locy)=value(locy)-gpi
      locy=lvn+nodplc(loc+19)
      value(locy)=value(locy)-gepr
      locy=lvn+nodplc(loc+20)
      value(locy)=value(locy)-go+geqcb
      locy=lvn+nodplc(loc+21)
      value(locy)=value(locy)-gpi-gm-geqcb
      locy=lvn+nodplc(loc+31)
      value(locy)=value(locy)+gccs
      locy=lvn+nodplc(loc+32)
      value(locy)=value(locy)-gccs
      locy=lvn+nodplc(loc+33)
      value(locy)=value(locy)-gccs
      locy=lvn+nodplc(loc+34)
      value(locy)=value(locy)-geqbx
      locy=lvn+nodplc(loc+35)
      value(locy)=value(locy)-geqbx
 1000 loc=nodplc(loc)
      go to 10
      end
