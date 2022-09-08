      subroutine load
      implicit double precision (a-h,o-z)
c
c     this routine zeroes-out and then loads the coefficient matrix.
c the active devices and the controlled sources are loaded by separate
c subroutines.
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
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
      dimension qcap(1),ccap(1)
      equivalence (qcap(1),value(1)),(ccap(1),value(2))
      dimension find(1),vind(1)
      equivalence (find(1),value(1)),(vind(1),value(2))
c
      call second(t1)
c
c  zero y matrix and current vector
c
      call zero8(value(lvn+1),nstop+nttbr)
c
c  resistors
c
      loc=locate(1)
   20 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 30
      locv=nodplc(loc+1)
      val=value(locv+1)
      locy=lvn+nodplc(loc+6)
      value(locy)=value(locy)+val
      locy=lvn+nodplc(loc+7)
      value(locy)=value(locy)+val
      locy=lvn+nodplc(loc+4)
      value(locy)=value(locy)-val
      locy=lvn+nodplc(loc+5)
      value(locy)=value(locy)-val
      loc=nodplc(loc)
      go to 20
c
c  capacitors
c
   30 loc=locate(2)
      if ((mode.eq.1).and.(modedc.ne.2)) go to 100
   40 if ((loc.eq.0).or.(nodplc(loc+12).ne.0)) go to 100
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      loct=nodplc(loc+8)
      ipoly=nodplc(loc+4)
      if (ipoly.eq.1) go to 43
      lcoef=nodplc(loc+7)
      call sizmem(nodplc(loc+7),ncoef)
   43 vcap=value(locv+2)
      if ((mode.eq.1).and.(initf.eq.2)) go to 45
      if ((nosolv.ne.0).and.(initf.eq.5)) go to 45
      vcap=value(lvnim1+node1)-value(lvnim1+node2)
   45 value(locv+3)=vcap
      if (mode.eq.1) go to 60
   47 if (initf.ne.6) go to 50
      qcap(lx0+loct)=qcap(lx1+loct)
      go to 60
   50 if (ipoly.eq.0) go to 53
      qcap(lx0+loct)=value(locv+1)*vcap
      if (initf.ne.5) go to 60
      if (nosolv.ne.0) qcap(lx0+loct)=value(locv+1)*value(locv+2)
      qcap(lx1+loct)=qcap(lx0+loct)
      go to 60
   53 call evpoly(qcap(lx0+loct),-1,lcoef,ncoef,locv+2,1,loc+8)
      if (initf.ne.5) go to 60
      if (nosolv.eq.0) go to 55
      vcap=value(locv+2)
      value(locv+3)=vcap
      call evpoly(qcap(lx0+loct),-1,lcoef,ncoef,locv+2,1,loc+8)
   55 qcap(lx1+loct)=qcap(lx0+loct)
   60 if (ipoly.eq.1) go to 62
      call evpoly(value(locv+1),0,lcoef,ncoef,locv+2,1,loc+8)
   62 if (mode.eq.1) go to 90
      call intgr8(geq,ceq,value(locv+1),loct)
      if (ipoly.eq.1) go to 65
      ceq=ceq-geq*vcap+ag(1)*qcap(lx0+loct)
   65 if(initf.ne.5) go to 70
      ccap(lx1+loct)=ccap(lx0+loct)
   70 locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)+geq
      locy=lvn+nodplc(loc+11)
      value(locy)=value(locy)+geq
      locy=lvn+nodplc(loc+5)
      value(locy)=value(locy)-geq
      locy=lvn+nodplc(loc+6)
      value(locy)=value(locy)-geq
      value(lvn+node1)=value(lvn+node1)-ceq
      value(lvn+node2)=value(lvn+node2)+ceq
   90 loc=nodplc(loc)
      go to 40
c
c  inductors
c
  100 if (jelcnt(3).eq.0) go to 400
      if (mode.eq.1) go to 150
      if (initf.eq.6) go to 150
      loc=locate(3)
  110 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 120
      locv=nodplc(loc+1)
      iptr=nodplc(loc+5)
      loct=nodplc(loc+11)
      ipoly=nodplc(loc+4)
      if (ipoly.eq.0) go to 115
      find(lx0+loct)=value(locv+1)*value(lvnim1+iptr)
      if ((initf.eq.5).and.(nosolv.ne.0))
     1   find(lx0+loct)=value(locv+1)*value(locv+2)
      go to 118
  115 lcoef=nodplc(loc+10)
      call sizmem(nodplc(loc+10),ncoef)
      cind=value(lvnim1+iptr)
      if ((initf.eq.5).and.(nosolv.ne.0)) cind=value(locv+2)
      value(locv+3)=cind
      call evpoly(find(lx0+loct),-1,lcoef,ncoef,locv+2,1,loc+11)
  118 loc=nodplc(loc)
      go to 110
  120 loc=locate(4)
  130 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 150
      locv=nodplc(loc+1)
      nl1=nodplc(loc+2)
      nl2=nodplc(loc+3)
      iptr1=nodplc(nl1+5)
      iptr2=nodplc(nl2+5)
      loct1=nodplc(nl1+11)
      loct2=nodplc(nl2+11)
      find(lx0+loct1)=find(lx0+loct1)+value(locv+1)*value(lvnim1+iptr2)
      find(lx0+loct2)=find(lx0+loct2)+value(locv+1)*value(lvnim1+iptr1)
      loc=nodplc(loc)
      go to 130
  150 loc=locate(3)
  160 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 300
      locv=nodplc(loc+1)
      iptr=nodplc(loc+5)
      loct=nodplc(loc+11)
      ipoly=nodplc(loc+4)
      if (ipoly.eq.1) go to 170
      lcoef=nodplc(loc+10)
      call sizmem(nodplc(loc+10),ncoef)
  170 cind=value(lvnim1+iptr)
      if ((nosolv.ne.0).and.(initf.eq.5)) cind=value(locv+2)
      value(locv+3)=cind
  180 if (mode.ne.1) go to 200
      veq=0.0d0
      req=0.0d0
      go to 210
  200 if (initf.ne.6) go to 205
      find(lx0+loct)=find(lx1+loct)
      go to 210
  205 if (initf.ne.5) go to 210
      find(lx1+loct)=find(lx0+loct)
  210 if (ipoly.eq.1) go to 220
      call evpoly(value(locv+1),0,lcoef,ncoef,locv+2,1,loc+11)
  220 if (mode.eq.1) go to 250
      call intgr8(req,veq,value(locv+1),loct)
      if (ipoly.eq.1) go to 250
      veq=veq-req*cind+ag(1)*find(lx0+loct)
  250 value(lvn+iptr)=veq
      if(initf.ne.5) go to 260
      vind(lx1+loct)=vind(lx0+loct)
  260 locy=lvn+nodplc(loc+13)
      value(locy)=-req
      locy=lvn+nodplc(loc+6)
      value(locy)=1.0d0
      locy=lvn+nodplc(loc+7)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+8)
      value(locy)=1.0d0
      locy=lvn+nodplc(loc+9)
      value(locy)=-1.0d0
      loc=nodplc(loc)
      go to 160
c
c  mutual inductances
c
  300 loc=locate(4)
  310 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 400
      locv=nodplc(loc+1)
      req=ag(1)*value(locv+1)
      locy=lvn+nodplc(loc+4)
      value(locy)=-req
      locy=lvn+nodplc(loc+5)
      value(locy)=-req
      loc=nodplc(loc)
      go to 310
c
c  nonlinear controlled sources
c
  400 call nlcsrc
c
c  voltage sources
c
      loc=locate(9)
  610 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 700
      locv=nodplc(loc+1)
      iptr=nodplc(loc+6)
      value(lvn+iptr)=value(locv+1)*sfactr
      locy=lvn+nodplc(loc+7)
      value(locy)=value(locy)+1.0d0
      locy=lvn+nodplc(loc+8)
      value(locy)=value(locy)-1.0d0
      locy=lvn+nodplc(loc+9)
      value(locy)=value(locy)+1.0d0
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)-1.0d0
      loc=nodplc(loc)
      go to 610
c
c  current sources
c
  700 loc=locate(10)
  710 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 800
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      val=value(locv+1)*sfactr
      value(lvn+node1)=value(lvn+node1)-val
      value(lvn+node2)=value(lvn+node2)+val
      loc=nodplc(loc)
      go to 710
c
c  call device model routines
c
  800 call diode
      call bjt
      call jfet
      call mosfet
c
c  transmission lines
c
      loc=locate(17)
  910 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 980
      locv=nodplc(loc+1)
      z0=value(locv+1)
      y0=1.0d0/z0
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      ibr1=nodplc(loc+8)
      ibr2=nodplc(loc+9)
      locy=lvn+nodplc(loc+10)
      value(locy)=value(locy)+y0
      locy=lvn+nodplc(loc+11)
      value(locy)=-y0
      locy=lvn+nodplc(loc+12)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+13)
      value(locy)=value(locy)+y0
      locy=lvn+nodplc(loc+14)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+15)
      value(locy)=-y0
      locy=lvn+nodplc(loc+16)
      value(locy)=+y0
      locy=lvn+nodplc(loc+17)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+18)
      value(locy)=+y0
      locy=lvn+nodplc(loc+19)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+20)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+23)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+27)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+28)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+31)
      value(locy)=-y0
      locy=lvn+nodplc(loc+32)
      value(locy)=-y0
      if (mode.ne.1) go to 920
      locy=lvn+nodplc(loc+21)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+22)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+24)
      value(locy)=-(1.0d0-gmin)*z0
      locy=lvn+nodplc(loc+25)
      value(locy)=-1.0d0
      locy=lvn+nodplc(loc+26)
      value(locy)=+1.0d0
      locy=lvn+nodplc(loc+29)
      value(locy)=-(1.0d0-gmin)*z0
      go to 950
  920 if (initf.ne.5) go to 930
      if (nosolv.ne.0) go to 925
      value(locv+3)=value(lvnim1+node3)-value(lvnim1+node4)
     1   +value(lvnim1+ibr2)*z0
      value(locv+4)=value(lvnim1+node1)-value(lvnim1+node2)
     1   +value(lvnim1+ibr1)*z0
      go to 930
  925 value(locv+3)=value(locv+7)+value(locv+8)*z0
      value(locv+4)=value(locv+5)+value(locv+6)*z0
  930 value(lvn+ibr1)=value(locv+3)
      value(lvn+ibr2)=value(locv+4)
  950 loc=nodplc(loc)
      go to 910
c
c  initialize nodes
c
  980 if(mode.ne.1) go to 995
      if(initf.ne.3.and.initf.ne.2) go to 995
      call sizmem(nsnod,nic)
      if(nic.eq.0) go to 995
      call sizmem(icnod,ntest)
      if(modedc.eq.2.and.ntest.ne.0) go to 995
      g=1.0d0
      do 990 i=1,nic
      locy=lvn+nodplc(nsmat+i)
      value(locy)=value(locy)+g
      node=nodplc(nsnod+i)
      value(lvn+node)=value(lvn+node)+value(nsval+i)*g
  990 continue
c
c  transient initial conditions (uic not specified)
c
  995 if(mode.ne.1) go to 1000
      if(modedc.ne.2) go to 1000
      if(nosolv.ne.0) go to 1000
      call sizmem(icnod,nic)
      if(nic.eq.0) go to 1000
      g=1.0d0
      do 996 i=1,nic
      locy=lvn+nodplc(icmat+i)
      value(locy)=value(locy)+g
      node=nodplc(icnod+i)
      value(lvn+node)=value(lvn+node)+value(icval+i)*g
  996 continue
c
c  finished
c
 1000 call second(t2)
      rstats(45)=rstats(45)+t2-t1
      return
      end
