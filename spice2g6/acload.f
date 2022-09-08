      subroutine acload
      implicit double precision (a-h,o-z)
c
c     this routine zeroes-out and then loads the complex coefficient
c     matrix
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
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      complex cval
c
c  zero y matrix and current vector
c
      call zero8(value(lvn+1),nstop+nttbr)
      call zero8(value(imvn+1),nstop+nttbr)
c
c  resistors
c
      loc=locate(1)
   20 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 30
      locv=nodplc(loc+1)
      val=value(locv+1)
      locy=lynl+nodplc(loc+6)
      value(locy)=value(locy)+val
      locy=lynl+nodplc(loc+7)
      value(locy)=value(locy)+val
      locy=lynl+nodplc(loc+4)
      value(locy)=value(locy)-val
      locy=lynl+nodplc(loc+5)
      value(locy)=value(locy)-val
      loc=nodplc(loc)
      go to 20
c
c  capacitors
c
   30 loc=locate(2)
   40 if ((loc.eq.0).or.(nodplc(loc+12).ne.0)) go to 50
      locv=nodplc(loc+1)
      val=omega*value(locv+1)
      locyi=imynl+nodplc(loc+10)
      value(locyi)=value(locyi)+val
      locyi=imynl+nodplc(loc+11)
      value(locyi)=value(locyi)+val
      locyi=imynl+nodplc(loc+5)
      value(locyi)=value(locyi)-val
      locyi=imynl+nodplc(loc+6)
      value(locyi)=value(locyi)-val
      loc=nodplc(loc)
      go to 40
c
c  inductors
c
   50 loc=locate(3)
   60 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 70
      locv=nodplc(loc+1)
      val=omega*value(locv+1)
      locyi=imynl+nodplc(loc+13)
      locy=lynl+nodplc(loc+13)
      value(locy)=0.0d0
      value(locyi)=-val
      locy=lynl+nodplc(loc+6)
      locyi=imynl+nodplc(loc+6)
      value(locy)=1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+7)
      locyi=imynl+nodplc(loc+7)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+8)
      locyi=imynl+nodplc(loc+8)
      value(locy)=1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+9)
      locyi=imynl+nodplc(loc+9)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      loc=nodplc(loc)
      go to 60
c
c  mutual inductors
c
   70 loc=locate(4)
   80 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 90
      locv=nodplc(loc+1)
      val=omega*value(locv+1)
      locy=lynl+nodplc(loc+4)
      locyi=imynl+nodplc(loc+4)
      value(locy)=0.0d0
      value(locyi)=-val
      locy=lynl+nodplc(loc+5)
      locyi=imynl+nodplc(loc+5)
      value(locy)=0.0d0
      value(locyi)=-val
      loc=nodplc(loc)
      go to 80
c
c  nonlinear voltage controlled current sources
c
   90 loc=locate(5)
   95 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 100
      ndim=nodplc(loc+4)
      lmat=nodplc(loc+7)
      loct=lx0+nodplc(loc+12)+2
      do 97 i=1,ndim
      val=value(loct)
      loct=loct+2
      locy=lynl+nodplc(lmat+1)
      value(locy)=value(locy)+val
      locy=lynl+nodplc(lmat+2)
      value(locy)=value(locy)-val
      locy=lynl+nodplc(lmat+3)
      value(locy)=value(locy)-val
      locy=lynl+nodplc(lmat+4)
      value(locy)=value(locy)+val
      lmat=lmat+4
   97 continue
      loc=nodplc(loc)
      go to 95
c
c  nonlinear voltage controlled voltage sources
c
  100 loc=locate(6)
  105 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 110
      ndim=nodplc(loc+4)
      lmat=nodplc(loc+8)
      loct=lx0+nodplc(loc+13)+3
      locy=lynl+nodplc(lmat+1)
      locyi=imynl+nodplc(lmat+1)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+2)
      locyi=imynl+nodplc(lmat+2)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+3)
      locyi=imynl+nodplc(lmat+3)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+4)
      locyi=imynl+nodplc(lmat+4)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      lmat=lmat+4
      do 107 i=1,ndim
      val=value(loct)
      loct=loct+2
      locy=lynl+nodplc(lmat+1)
      value(locy)=value(locy)-val
      locy=lynl+nodplc(lmat+2)
      value(locy)=value(locy)+val
      lmat=lmat+2
  107 continue
      loc=nodplc(loc)
      go to 105
c
c  nonlinear current controlled current sources
c
  110 loc=locate(7)
  115 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 120
      ndim=nodplc(loc+4)
      lmat=nodplc(loc+7)
      loct=lx0+nodplc(loc+12)+2
      do 117 i=1,ndim
      val=value(loct)
      loct=loct+2
      locy=lynl+nodplc(lmat+1)
      locyi=imynl+nodplc(lmat+1)
      value(locy)=+val
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+2)
      locyi=imynl+nodplc(lmat+2)
      value(locy)=-val
      value(locyi)=0.0d0
      lmat=lmat+2
  117 continue
      loc=nodplc(loc)
      go to 115
c
c  nonlinear current controlled voltage sources
c
  120 loc=locate(8)
  125 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 140
      ndim=nodplc(loc+4)
      lmat=nodplc(loc+8)
      loct=lx0+nodplc(loc+13)+3
      locy=lynl+nodplc(lmat+1)
      locyi=imynl+nodplc(lmat+1)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+2)
      locyi=imynl+nodplc(lmat+2)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+3)
      locyi=imynl+nodplc(lmat+3)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(lmat+4)
      locyi=imynl+nodplc(lmat+4)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      lmat=lmat+4
      do 127 i=1,ndim
      val=value(loct)
      loct=loct+2
      locy=lynl+nodplc(lmat+i)
      value(locy)=value(locy)-val
  127 continue
      loc=nodplc(loc)
      go to 125
c
c  voltage sources
c
  140 loc=locate(9)
  150 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 160
      locv=nodplc(loc+1)
      iptr=nodplc(loc+6)
      value(lvn+iptr)=value(locv+2)
      value(imvn+iptr)=value(locv+3)
      locy=lynl+nodplc(loc+7)
      value(locy)=value(locy)+1.0d0
      locy=lynl+nodplc(loc+8)
      value(locy)=value(locy)-1.0d0
      locy=lynl+nodplc(loc+9)
      value(locy)=value(locy)+1.0d0
      locy=lynl+nodplc(loc+10)
      value(locy)=value(locy)-1.0d0
      loc=nodplc(loc)
      go to 150
c
c  current sources
c
  160 loc=locate(10)
  170 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 200
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      value(lvn+node1)=value(lvn+node1)-value(locv+2)
      value(imvn+node1)=value(imvn+node1)-value(locv+3)
      value(lvn+node2)=value(lvn+node2)+value(locv+2)
      value(imvn+node2)=value(imvn+node2)+value(locv+3)
      loc=nodplc(loc)
      go to 170
c
c  diodes
c
  200 loc=locate(11)
  210 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 250
      locv=nodplc(loc+1)
      area=value(locv+1)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+11)
      gspr=value(locm+2)*area
      geq=value(loct+2)
      xceq=value(loct+4)*omega
      locy=lynl+nodplc(loc+13)
      value(locy)=value(locy)+gspr
      locy=lynl+nodplc(loc+14)
      locyi=imynl+nodplc(loc+14)
      value(locy)=value(locy)+geq
      value(locyi)=value(locyi)+xceq
      locy=lynl+nodplc(loc+15)
      locyi=imynl+nodplc(loc+15)
      value(locy)=value(locy)+geq+gspr
      value(locyi)=value(locyi)+xceq
      locy=lynl+nodplc(loc+7)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+8)
      locyi=imynl+nodplc(loc+8)
      value(locy)=value(locy)-geq
      value(locyi)=value(locyi)-xceq
      locy=lynl+nodplc(loc+9)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+10)
      locyi=imynl+nodplc(loc+10)
      value(locy)=value(locy)-geq
      value(locyi)=value(locyi)-xceq
      loc=nodplc(loc)
      go to 210
c
c  bjts
c
  250 loc=locate(12)
  260 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 300
      locv=nodplc(loc+1)
      area=value(locv+1)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+22)
      gcpr=value(locm+20)*area
      gepr=value(locm+19)*area
      gpi=value(loct+4)
      gmu=value(loct+5)
      gm=value(loct+6)
      go=value(loct+7)
      xgm=0.0d0
      td=value(locm+28)
      if(td.eq.0.0d0) go to 270
      arg=td*omega
      gm=gm+go
      xgm=-gm*dsin(arg)
      gm=gm*dcos(arg)-go
  270 gx=value(loct+16)
      xcpi=value(loct+9)*omega
      xcmu=value(loct+11)*omega
      xcbx=value(loct+15)*omega
      xccs=value(loct+13)*omega
      xcmcb=value(loct+17)*omega
      locy=lynl+nodplc(loc+24)
      value(locy)=value(locy)+gcpr
      locy=lynl+nodplc(loc+25)
      locyi=imynl+nodplc(loc+25)
      value(locy)=value(locy)+gx
      value(locyi)=value(locyi)+xcbx
      locy=lynl+nodplc(loc+26)
      value(locy)=value(locy)+gepr
      locy=lynl+nodplc(loc+27)
      locyi=imynl+nodplc(loc+27)
      value(locy)=value(locy)+gmu+go+gcpr
      value(locyi)=value(locyi)+xcmu+xccs+xcbx
      locy=lynl+nodplc(loc+28)
      locyi=imynl+nodplc(loc+28)
      value(locy)=value(locy)+gx+gpi+gmu
      value(locyi)=value(locyi)+xcpi+xcmu+xcmcb
      locy=lynl+nodplc(loc+29)
      locyi=imynl+nodplc(loc+29)
      value(locy)=value(locy)+gpi+gepr+gm+go
      value(locyi)=value(locyi)+xcpi+xgm
      locy=lynl+nodplc(loc+10)
      value(locy)=value(locy)-gcpr
      locy=lynl+nodplc(loc+11)
      value(locy)=value(locy)-gx
      locy=lynl+nodplc(loc+12)
      value(locy)=value(locy)-gepr
      locy=lynl+nodplc(loc+13)
      value(locy)=value(locy)-gcpr
      locy=lynl+nodplc(loc+14)
      locyi=imynl+nodplc(loc+14)
      value(locy)=value(locy)-gmu+gm
      value(locyi)=value(locyi)-xcmu+xgm
      locy=lynl+nodplc(loc+15)
      locyi=imynl+nodplc(loc+15)
      value(locy)=value(locy)-gm-go
      value(locyi)=value(locyi)-xgm
      locy=lynl+nodplc(loc+16)
      value(locy)=value(locy)-gx
      locy=lynl+nodplc(loc+17)
      locyi=imynl+nodplc(loc+17)
      value(locy)=value(locy)-gmu
      value(locyi)=value(locyi)-xcmu-xcmcb
      locy=lynl+nodplc(loc+18)
      locyi=imynl+nodplc(loc+18)
      value(locy)=value(locy)-gpi
      value(locyi)=value(locyi)-xcpi
      locy=lynl+nodplc(loc+19)
      value(locy)=value(locy)-gepr
      locy=lynl+nodplc(loc+20)
      locyi=imynl+nodplc(loc+20)
      value(locy)=value(locy)-go
      value(locyi)=value(locyi)+xcmcb
      locy=lynl+nodplc(loc+21)
      locyi=imynl+nodplc(loc+21)
      value(locy)=value(locy)-gpi-gm
      value(locyi)=value(locyi)-xcpi-xgm-xcmcb
      locyi=imynl+nodplc(loc+31)
      value(locyi)=value(locyi)+xccs
      locyi=imynl+nodplc(loc+32)
      value(locyi)=value(locyi)-xccs
      locyi=imynl+nodplc(loc+33)
      value(locyi)=value(locyi)-xccs
      locyi=imynl+nodplc(loc+34)
      value(locyi)=value(locyi)-xcbx
      locyi=imynl+nodplc(loc+35)
      value(locyi)=value(locyi)-xcbx
      loc=nodplc(loc)
      go to 260
c
c  jfets
c
  300 loc=locate(13)
  310 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 350
      locv=nodplc(loc+1)
      area=value(locv+1)
      locm=nodplc(loc+7)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+19)
      gdpr=value(locm+4)*area
      gspr=value(locm+5)*area
      gm=value(loct+5)
      gds=value(loct+6)
      ggs=value(loct+7)
      xgs=value(loct+9)*omega
      ggd=value(loct+8)
      xgd=value(loct+11)*omega
      locy=lynl+nodplc(loc+20)
      value(locy)=value(locy)+gdpr
      locy=lynl+nodplc(loc+21)
      locyi=imynl+nodplc(loc+21)
      value(locy)=value(locy)+ggd+ggs
      value(locyi)=value(locyi)+xgd+xgs
      locy=lynl+nodplc(loc+22)
      value(locy)=value(locy)+gspr
      locy=lynl+nodplc(loc+23)
      locyi=imynl+nodplc(loc+23)
      value(locy)=value(locy)+gdpr+gds+ggd
      value(locyi)=value(locyi)+xgd
      locy=lynl+nodplc(loc+24)
      locyi=imynl+nodplc(loc+24)
      value(locy)=value(locy)+gspr+gds+gm+ggs
      value(locyi)=value(locyi)+xgs
      locy=lynl+nodplc(loc+9)
      value(locy)=value(locy)-gdpr
      locy=lynl+nodplc(loc+10)
      locyi=imynl+nodplc(loc+10)
      value(locy)=value(locy)-ggd
      value(locyi)=value(locyi)-xgd
      locy=lynl+nodplc(loc+11)
      locyi=imynl+nodplc(loc+11)
      value(locy)=value(locy)-ggs
      value(locyi)=value(locyi)-xgs
      locy=lynl+nodplc(loc+12)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+13)
      value(locy)=value(locy)-gdpr
      locy=lynl+nodplc(loc+14)
      locyi=imynl+nodplc(loc+14)
      value(locy)=value(locy)-ggd+gm
      value(locyi)=value(locyi)-xgd
      locy=lynl+nodplc(loc+15)
      value(locy)=value(locy)-gds-gm
      locy=lynl+nodplc(loc+16)
      locyi=imynl+nodplc(loc+16)
      value(locy)=value(locy)-ggs-gm
      value(locyi)=value(locyi)-xgs
      locy=lynl+nodplc(loc+17)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+18)
      value(locy)=value(locy)-gds
      loc=nodplc(loc)
      go to 310
c
c  mosfets
c
  350 loc=locate(14)
  360 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 400
      locv=nodplc(loc+1)
      locm=nodplc(loc+8)
      itype=nodplc(locm+2)
      locm=nodplc(locm+1)
      devmod=value(locv+8)
      xnrm=1.0d0
      xrev=0.0d0
      if (devmod.ge.0.0d0) go to 370
      xnrm=0.0d0
      xrev=1.0d0
  370 loct=lx0+nodplc(loc+26)
      if (value(locm+7).eq.0.0d0.and.
     1   value(locm+8).eq.0.0d0) go to 375
      gdpr=value(locm+7)
      gspr=value(locm+8)
      go to 380
  375 gdpr=value(locm+16)/value(locv+13)
      gspr=value(locm+16)/value(locv+14)
  380 gm=value(loct+7)
      gds=value(loct+8)
      gmbs=value(loct+9)
      gbd=value(loct+10)
      gbs=value(loct+11)
      capbd=value(loct+24)
      capbs=value(loct+26)
cc
cc    charge oriented model parameters
cc
      xl=value(locv+1)-2.0d0*value(locm+28)
      xw=value(locv+2)
      xqco=value(locm+35)
      xqc=value(locv+15)
      covlgs=value(locm+13)*xw
      covlgd=value(locm+14)*xw
      covlgb=value(locm+15)*xl
      if (xqco.gt.0.5d0) go to 385
      cggb=value(loct+18)
      cgdb=value(loct+19)
      cgsb=value(loct+20)
      cbgb=value(loct+21)
      cbdb=value(loct+22)
      cbsb=value(loct+23)
      gcg=(cggb+cbgb)*omega
      gcd=(cgdb+cbdb)*omega
      gcs=(cgsb+cbsb)*omega
      xcgxd=-xqc*gcg
      xcgxs=-(1.0d0-xqc)*gcg
      xcdxd=-xqc*gcd
      xcdxs=-(1.0d0-xqc)*gcd
      xcsxd=-xqc*gcs
      xcsxs=-(1.0d0-xqc)*gcs
      xcdgb=xcgxd-covlgd*omega
      xcddb=xcdxd+(capbd+covlgd)*omega
      xcdsb=xcsxd
      xcsgb=xcgxs-covlgs*omega
      xcsdb=xcdxs
      xcssb=xcsxs+(capbs+covlgs)*omega
      xcggb=(cggb+covlgd+covlgs+covlgb)*omega
      xcgdb=(cgdb-covlgd)*omega
      xcgsb=(cgsb-covlgs)*omega
      xcbgb=(cbgb-covlgb)*omega
      xcbdb=(cbdb-capbd)*omega
      xcbsb=(cbsb-capbs)*omega
      go to 390
c
c     meyer"s model parameters
c
  385 xcgs=(value(loct+12)+covlgs)*omega
      xcgd=(value(loct+14)+covlgd)*omega
      xcgb=(value(loct+16)+covlgb)*omega
      xbd=capbd*omega
      xbs=capbs*omega
cc
cc    do the mapping from meyer"s model into charge oriented model
cc
      xcggb=xcgd+xcgs+xcgb
      xcgdb=-xcgd
      xcgsb=-xcgs
      xcbgb=-xcgb
      xcbdb=-xbd
      xcbsb=-xbs
      xcddb=xcgd+xbd
      xcssb=xcgs+xbs
cc    xcgsb=-xcgb
      xcdgb=-xcgd
      xcsgb=-xcgs
      xcdsb=0.0d0
      xcsdb=0.0d0
cc
  390 locyi=imynl+nodplc(loc+28)
      value(locyi)=value(locyi)+xcggb
      locyi=imynl+nodplc(loc+30)
      value(locyi)=value(locyi)-xcbgb-xcbdb-xcbsb
      locyi=imynl+nodplc(loc+31)
      value(locyi)=value(locyi)+xcddb
      locyi=imynl+nodplc(loc+32)
      value(locyi)=value(locyi)+xcssb
      locyi=imynl+nodplc(loc+11)
      value(locyi)=value(locyi)-xcggb-xcgdb-xcgsb
      locyi=imynl+nodplc(loc+12)
      value(locyi)=value(locyi)+xcgdb
      locyi=imynl+nodplc(loc+13)
      value(locyi)=value(locyi)+xcgsb
      locyi=imynl+nodplc(loc+15)
      value(locyi)=value(locyi)+xcbgb
      locyi=imynl+nodplc(loc+16)
      value(locyi)=value(locyi)+xcbdb
      locyi=imynl+nodplc(loc+17)
      value(locyi)=value(locyi)+xcbsb
      locyi=imynl+nodplc(loc+19)
      value(locyi)=value(locyi)+xcdgb
      locyi=imynl+nodplc(loc+20)
      value(locyi)=value(locyi)-xcdgb-xcddb-xcdsb
      locyi=imynl+nodplc(loc+21)
      value(locyi)=value(locyi)+xcdsb
      locyi=imynl+nodplc(loc+22)
      value(locyi)=value(locyi)+xcsgb
      locyi=imynl+nodplc(loc+24)
      value(locyi)=value(locyi)-xcsgb-xcsdb-xcssb
      locyi=imynl+nodplc(loc+25)
      value(locyi)=value(locyi)+xcsdb
      locy=lynl+nodplc(loc+27)
      value(locy)=value(locy)+gdpr
      locy=lynl+nodplc(loc+29)
      value(locy)=value(locy)+gspr
      locy=lynl+nodplc(loc+30)
      value(locy)=value(locy)+gbd+gbs
      locy=lynl+nodplc(loc+31)
      value(locy)=value(locy)+gdpr+gds+gbd+xrev*(gm+gmbs)
      locy=lynl+nodplc(loc+32)
      value(locy)=value(locy)+gspr+gds+gbs+xnrm*(gm+gmbs)
      locy=lynl+nodplc(loc+10)
      value(locy)=value(locy)-gdpr
      locy=lynl+nodplc(loc+14)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+16)
      value(locy)=value(locy)-gbd
      locy=lynl+nodplc(loc+17)
      value(locy)=value(locy)-gbs
      locy=lynl+nodplc(loc+18)
      value(locy)=value(locy)-gdpr
      locy=lynl+nodplc(loc+19)
      value(locy)=value(locy)+(xnrm-xrev)*gm
      locy=lynl+nodplc(loc+20)
      value(locy)=value(locy)-gbd+(xnrm-xrev)*gmbs
      locy=lynl+nodplc(loc+21)
      value(locy)=value(locy)-gds-xnrm*(gm+gmbs)
      locy=lynl+nodplc(loc+22)
      value(locy)=value(locy)-(xnrm-xrev)*gm
      locy=lynl+nodplc(loc+23)
      value(locy)=value(locy)-gspr
      locy=lynl+nodplc(loc+24)
      value(locy)=value(locy)-gbs-(xnrm-xrev)*gmbs
      locy=lynl+nodplc(loc+25)
      value(locy)=value(locy)-gds-xrev*(gm+gmbs)
      loc=nodplc(loc)
      go to 360
c
c  transmission lines
c
  400 loc=locate(17)
  410 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 1000
      locv=nodplc(loc+1)
      z0=value(locv+1)
      y0=1.0d0/z0
      td=value(locv+2)
      arg=-omega*td
      rval=dcos(arg)
      xval=dsin(arg)
      locy=lynl+nodplc(loc+10)
      value(locy)=value(locy)+y0
      locy=lynl+nodplc(loc+11)
      locyi=imynl+nodplc(loc+11)
      value(locy)=-y0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+12)
      locyi=imynl+nodplc(loc+12)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+13)
      value(locy)=value(locy)+y0
      locy=lynl+nodplc(loc+14)
      locyi=imynl+nodplc(loc+14)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+15)
      locyi=imynl+nodplc(loc+15)
      value(locy)=-y0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+16)
      locyi=imynl+nodplc(loc+16)
      value(locy)=+y0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+17)
      locyi=imynl+nodplc(loc+17)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+18)
      locyi=imynl+nodplc(loc+18)
      value(locy)=+y0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+19)
      locyi=imynl+nodplc(loc+19)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+20)
      locyi=imynl+nodplc(loc+20)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+21)
      locyi=imynl+nodplc(loc+21)
      value(locy)=-rval
      value(locyi)=-xval
      locy=lynl+nodplc(loc+22)
      locyi=imynl+nodplc(loc+22)
      value(locy)=+rval
      value(locyi)=+xval
      locy=lynl+nodplc(loc+23)
      locyi=imynl+nodplc(loc+23)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+24)
      locyi=imynl+nodplc(loc+24)
      value(locy)=-rval*z0
      value(locyi)=-xval*z0
      locy=lynl+nodplc(loc+25)
      locyi=imynl+nodplc(loc+25)
      value(locy)=-rval
      value(locyi)=-xval
      locy=lynl+nodplc(loc+26)
      locyi=imynl+nodplc(loc+26)
      value(locy)=+rval
      value(locyi)=+xval
      locy=lynl+nodplc(loc+27)
      locyi=imynl+nodplc(loc+27)
      value(locy)=-1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+28)
      locyi=imynl+nodplc(loc+28)
      value(locy)=+1.0d0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+29)
      locyi=imynl+nodplc(loc+29)
      value(locy)=-rval*z0
      value(locyi)=-xval*z0
      locy=lynl+nodplc(loc+31)
      locyi=imynl+nodplc(loc+31)
      value(locy)=-y0
      value(locyi)=0.0d0
      locy=lynl+nodplc(loc+32)
      locyi=imynl+nodplc(loc+32)
      value(locy)=-y0
      value(locyi)=0.0d0
      loc=nodplc(loc)
      go to 410
c
c  finished
c
 1000 return
      end
