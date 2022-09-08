      subroutine dinit
      implicit double precision (a-h,o-z)
c
c     this routine performs storage-allocation and one-time computation
c needed to do the small-signal distortion analysis.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      call getm8(ld0,ndist)
      call getm16(ld1,5*nstop)
c
c  bipolar junction transistors
c
      loc=locate(12)
  100 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 200
      locv=nodplc(loc+1)
      area=value(locv+1)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+22)
      locd=ld0+nodplc(loc+23)
      csat=value(locm+1)*area
      ova=value(locm+4)
      tf=value(locm+24)
      tr=value(locm+33)
      czbe=value(locm+21)*area
      czbc=value(locm+29)*area
      pe=value(locm+22)
      xme=value(locm+23)
      pc=value(locm+30)
      xmc=value(locm+31)
      fcpe=value(locm+46)
      fcpc=value(locm+50)
      vbe=value(loct)
      vbc=value(loct+1)
      gpi=value(loct+4)
      go=value(loct+7)
      gm=value(loct+6)
      gmu=value(loct+5)
      if (vbe.gt.0.0d0) go to 110
      evbe=1.0d0
      cbe=csat*vbe/vt
      go to 120
  110 evbe=dexp(vbe/vt)
      cbe=csat*(evbe-1.0d0)
  120 if (vbc.gt.0.0d0) go to 130
      evbc=1.0d0
      cbc=csat*vbc/vt
      arg=1.0d0-vbc/pc
      go to 140
  130 evbc=dexp(vbc/vt)
      cbc=csat*(evbc-1.0d0)
  140 if (vbe.ge.fcpe) go to 150
      arg=1.0d0-vbe/pe
      sarg=dexp(xme*dlog(arg))
      cjeo=czbe/sarg
      argbe=pe-vbe
      cje1=xme*cjeo/argbe
      cje2=(1.0d0+xme)*cje1/argbe
      go to 160
  150 denom=dexp((1.0d0+xme)*dlog(1.0d0-fcpe))
      cjeo=czbe*(1.0d0-fcpe*(1.0d0+xme)+xme*vbe/pe)/denom
      cje1=czbe*xme/(denom*pe)
      cje2=0.0d0
  160 if (vbc.ge.fcpc) go to 170
      arg=1.0d0-vbc/pc
      sarg=dexp(xmc*dlog(arg))
      cjco=czbc/sarg
      argbc=pc-vbc
      cjc1=xmc*cjco/argbc
      cjc2=(1.0d0+xmc)*cjc1/argbc
      go to 180
  170 denom=dexp((1.0d0+xmc)*dlog(1.0d0-fcpc))
      cjco=czbc*(1.0d0-fcpc*(1.0d0+xmc)+xmc*vbc/pc)/denom
      cjc1=czbc*xmc/(denom*pc)
      cjc2=0.0d0
  180 twovt=vt+vt
      go2=(-go+csat*(evbe+evbc)*ova)/twovt
      gmo2=(cbe+csat)*ova/vt-2.0d0*go2
      gm2=(gm+go)/twovt-gmo2-go2
      gmu2=gmu/twovt
      if (vbc.le.0.0d0) gmu2=0.0d0
      gpi2=gpi/twovt
      if (vbe.le.0.0d0) gpi2=0.0d0
      cbo=tf*csat*evbe/vt
      cbor=tr*csat*evbc/vt
      cb1=cbo/vt
      cb1r=cbor/vt
      trivt=3.0d0*vt
      go3=-(go2+(cbc+csat)*ova/twovt)/trivt
      gmo23=-3.0d0*go3
      gm2o3=-gmo23+(cbe+csat)*ova/(vt*twovt)
      gm3=(gm2-(cbe-cbc)*ova/twovt)/trivt
      gmu3=gmu2/trivt
      gpi3=gpi2/trivt
      cb2=cb1/twovt
      cb2r=cb1r/twovt
      value(locd)=cje1
      value(locd+1)=cje2
      value(locd+2)=cjc1
      value(locd+3)=cjc2
      value(locd+4)=go2
      value(locd+5)=gmo2
      value(locd+6)=gm2
      value(locd+7)=gmu2
      value(locd+8)=gpi2
      value(locd+9)=cbo
      value(locd+10)=cbor
      value(locd+11)=cb1
      value(locd+12)=cb1r
      value(locd+13)=go3
      value(locd+14)=gmo23
      value(locd+15)=gm2o3
      value(locd+16)=gm3
      value(locd+17)=gmu3
      value(locd+18)=gpi3
      value(locd+19)=cb2
      value(locd+20)=cb2r
      loc=nodplc(loc)
      go to 100
c
c  diodes
c
  200 loc=locate(11)
  210 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 300
      locv=nodplc(loc+1)
      area=value(locv+1)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+11)
      locd=ld0+nodplc(loc+12)
      csat=value(locm+1)*area
      vte=value(locm+3)*vt
      tau=value(locm+4)
      czero=value(locm+5)*area
      phib=value(locm+6)
      xm=value(locm+7)
      fcpb=value(locm+12)
      vd=value(loct)
      geq=value(loct+2)
      evd=1.0d0
      if (vd.ge.0.0d0) evd=dexp(vd/vte)
      if (vd.ge.fcpb) go to 220
      arg=1.0d0-vd/phib
      sarg=dexp(xm*dlog(arg))
      cdjo=czero/sarg
      argd=phib-vd
      cdj1=xm*cdjo/argd
      cdj2=(1.0d0+xm)*cdj1/argd
      go to 230
  220 denom=dexp((1.0d0+xm)*dlog(1.0d0-fcpb))
      cdjo=czero*(1.0d0-fcpb*(1.0d0+xm)+xm*vd/phib)/denom
      cdj1=czero*xm/(denom*phib)
      cdj2=0.0d0
      cdj2=0.0d0
  230 cdbo=tau*csat*evd/vte
      cdb1=cdbo/vte
      twovte=2.0d0*vte
      geq2=geq/twovte
      if (vd.le.0.0d0) geq2=0.0d0
      trivte=3.0d0*vte
      geq3=geq2/trivte
      cdb2=cdb1/twovte
      value(locd)=cdj1
      value(locd+1)=cdj2
      value(locd+2)=cdbo
      value(locd+3)=cdb1
      value(locd+4)=geq2
      value(locd+5)=geq3
      value(locd+6)=cdb2
      loc=nodplc(loc)
      go to 210
c
c  finished
c
  300 return
      end
