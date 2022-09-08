      subroutine modchk
      implicit double precision (a-h,o-z)
c
c     this routine performs one-time processing of device model para-
c meters and prints out a device model summary.  it also reserves the
c additional nodes required by nonzero device extrinsic resistances.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension itab(50),atable(12)
      dimension cpar(2),btable(12)
      dimension antype(4),aptype(4)
      dimension ipar(5),ampar(115),defval(115),ifmt(115),ivchk(115)
      dimension titled(4),titleb(4),titlej(4),titlem(4)
      data titled / 8hdiode mo, 8hdel para, 8hmeters  , 8h         /
      data titleb / 8hbjt mode, 8hl parame, 8hters    , 8h         /
      data titlej / 8hjfet mod, 8hel param, 8heters   , 8h         /
      data titlem / 8hmosfet m, 8hodel par, 8hameters , 8h         /
      data antype /1h ,3hnpn,3hnjf,4hnmos/
      data aptype /1h ,3hpnp,3hpjf,4hpmos/
      data ipar /0,14,60,72,114/
      data cpar / 3hc2 ,3hc4 /
      data aundef /2h.u/
      data ampar /
     1   6his    ,6hrs    ,6hn     ,6htt    ,6hcjo   ,6hvj    ,6hm     ,
     2   6heg    ,6hxti   ,6hkf    ,6haf    ,6hfc    ,6hbv    ,6hibv   ,
     1   6his    ,6hbf    ,6hnf    ,6hvaf   ,6hikf   ,6hise   ,6hne    ,
     2   6hbr    ,6hnr    ,6hvar   ,6hikr   ,6hisc   ,6hnc    ,6h0     ,
     3   6h0     ,6hrb    ,6hirb   ,6hrbm   ,6hre    ,6hrc    ,6hcje   ,
     4   6hvje   ,6hmje   ,6htf    ,6hxtf   ,6hvtf   ,6hitf   ,6hptf   ,
     5   6hcjc   ,6hvjc   ,6hmjc   ,6hxcjc  ,6htr    ,6h0     ,6h0     ,
     6   6h0     ,6h0     ,6hcjs   ,6hvjs   ,6hmjs   ,6hxtb   ,6heg    ,
     7   6hxti   ,6hkf    ,6haf    ,6hfc    ,
     1   6hvto   ,6hbeta  ,6hlambda,6hrd    ,6hrs    ,6hcgs   ,6hcgd   ,
     2   6hpb    ,6his    ,6hkf    ,6haf    ,6hfc    ,
     1   6hlevel ,6hvto   ,6hkp    ,6hgamma ,6hphi   ,6hlambda,6hrd    ,
     2   6hrs    ,6hcbd   ,6hcbs   ,6his    ,6hpb    ,6hcgso  ,6hcgdo  ,
     3   6hcgbo  ,6hrsh   ,6hcj    ,6hmj    ,6hcjsw  ,6hmjsw  ,6hjs    ,
     4   6htox   ,6hnsub  ,6hnss   ,6hnfs   ,6htpg   ,6hxj    ,6hld    ,
     5   6huo    ,6hucrit ,6huexp  ,6hutra  ,6hvmax  ,6hneff  ,6hxqc   ,
     6   6hkf    ,6haf    ,6hfc    ,6hdelta ,6htheta ,6heta   ,6hkappa ,
     7   0.0d0   /
      data defval /
     1   1.0d-14,  0.0d0,  1.0d0,2*0.0d0,  1.0d0,  0.5d0, 1.11d0,
     2     3.0d0,  0.0d0,  1.0d0,  0.5d0,  0.0d0, 1.0d-3,
     1   1.0d-16,100.0d0,  1.0d0,3*0.0d0,  1.5d0,2*1.0d0,3*0.0d0,
     2     2.0d0,  0.0d0,  1.0d0,6*0.0d0, 0.75d0, 0.33d0,2*0.0d0,
     3   4*0.0d0, 0.75d0, 0.33d0,  1.0d0,6*0.0d0, 0.75d0,2*0.0d0,
     4    1.11d0,  3.0d0,  0.0d0,  1.0d0,  0.5d0,
     1    -2.0d0, 1.0d-4,5*0.0d0,  1.0d0,1.0d-14,  0.0d0,  1.0d0,
     2     0.5d0,
     1     1.0d0,  0.0d0, 2.0d-5,  0.0d0,  0.6d0,5*0.0d0,1.0d-14,
     2     0.8d0,5*0.0d0,  0.5d0,  0.0d0, 0.33d0,5*0.0d0,  1.0d0,
     3   2*0.0d0,600.0d0, 1.0d+4,3*0.0d0,  1.0d0,  1.0d0,  0.0d0,
     4     1.0d0,  0.5d0,3*0.0d0,  0.2d0,
     5     0.0d0/
      data ifmt /
     1   4,1,1,2,2,1,1,1,1,2,1,1,2,2,
     2   4,3,3,2,2,2,1,3,3,2,2,2,1,0,0,1,2,1,1,1,2,1,1,2,2,2,2,1,2,1,
     2   1,1,2,0,0,0,0,2,1,1,2,1,1,2,2,2,
     3   3,4,1,1,1,2,2,1,2,2,1,1,
     4   3,3,4,1,1,2,1,1,2,2,2,1,2,2,2,1,2,1,2,1,2,2,2,2,2,1,2,2,
     4   1,2,1,1,2,1,1,2,1,1,1,1,1,1,
     5   0/
      data ivchk /
     1   0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     2   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     2   0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,
     3   -1,0,0,0,0,0,0,0,0,0,0,0,
     4   0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,-1,
     4   0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     5   0/
c
c
      tnom=value(itemps+1)+ctok
      xkt=boltz*tnom
      vt=xkt/charge
      xni=1.45d16
      egfet=1.16d0-(7.02d-4*tnom*tnom)/(tnom+1108.0d0)
      nummod=jelcnt(21)+jelcnt(22)+jelcnt(23)+jelcnt(24)
      if (nummod.eq.0) go to 1000
c
c  special preprocessing for mosfet models
c
      loc=locate(24)
    5 if (loc.eq.0) go to 35
      locv=nodplc(loc+1)
      type=nodplc(loc+2)
c
c     default parameters for higher level mos models
c
      lev=value(locv+1)
      if (value(locv+1).eq.aundef) lev=1
      if (value(locv+23).ne.aundef) xnsub=value(locv+23)*1.0d6
      if (value(locv+22).eq.aundef.and.lev.gt.1) value(locv+22)=1.0d-7
      if (value(locv+22).eq.aundef) go to 33
      cox=epsox/value(locv+22)
c
c     compute kp, if not input, using default mobility 600 cm**2/v*sec
c
      if (value(locv+3).ne.aundef) go to 10
      if (value(locv+29).eq.aundef) value(locv+29)=600.0d0
      value(locv+3)=value(locv+29)*cox*1.0d-4
   10 if (value(locv+23).eq.aundef) go to 33
      if (xnsub.le.xni) go to 30
c
c     nsub nonzero => process oriented model
c
      if (value(locv+5).eq.aundef) value(locv+5)=
     1   dmax1((2.0d0*vt*dlog(xnsub/xni)),0.1d0)
      fermis=type*0.5d0*value(locv+5)
      wkfng=3.2d0
      if (value(locv+26).eq.aundef) value(locv+26)=1.0d0
      if (value(locv+26).eq.0.0d0) go to 15
c
c  polysilicon gate
c
      fermig=type*value(locv+26)*0.5d0*egfet
      wkfng=3.25d0+0.5d0*egfet-fermig
   15 wkfngs=wkfng-(3.25d0+0.5d0*egfet+fermis)
      if (value(locv+4).eq.aundef)
     1   value(locv+4)=dsqrt(2.0d0*epssil*charge*xnsub)/cox
c
c     computed vto
c
      if (value(locv+2).ne.aundef) go to 20
      if (value(locv+24).eq.aundef) value(locv+24)=0.0d0
      value(locv+44)=wkfngs-value(locv+24)*1.0d4*charge/cox
      value(locv+2)=value(locv+44)
     1   +type*(value(locv+4)*dsqrt(value(locv+5))+value(locv+5))
      go to 25
c
c     measured vto has been input
c
   20 value(locv+44)=value(locv+2)
     1   -type*(value(locv+4)*dsqrt(value(locv+5))+value(locv+5))
   25 value(locv+45)=dsqrt((epssil+epssil)/(charge*xnsub))
      go to 33
   30 value(locv+23)=0.0d0
      write (iofile,31) value(locv)
   31 format('0*error*:  nsub <= ni in mosfet model ',a8,/)
      nogo=1
c
c   special processing for mos3: limit kappa>0,
c   set to zero lambda,ucrit,uexp and utar
c
   33 if (lev.ne.3) go to 34
      if (value(locv+42).eq.aundef) value(locv+42)=0.2d0
      value(locv+6)=0.0d0
      value(locv+30)=0.0d0
      value(locv+31)=0.0d0
      value(locv+32)=0.0d0
   34 loc=nodplc(loc)
      go to 5
c
c     cycle thru devices
c
   35 kntlim=lwidth/11
      do 390 id=1,4
      if (jelcnt(id+20).eq.0) go to 390
      locm=ipar(id)
      nopar=ipar(id+1)-locm
      do 45 i=1,nopar
      if (ifmt(locm+i).ge.3) go to 40
      itab(i)=0
      go to 45
   40 itab(i)=ifmt(locm+i)-2
   45 continue
c
c  assign default values
c
      loc=locate(id+20)
   50 if (loc.eq.0) go to 70
      locv=nodplc(loc+1)
      do 65 i=1,nopar
      if (value(locv+i).eq.aundef) go to 62
      if (ivchk(locm+i).lt.0) go to 55
      if (value(locv+i).lt.0.0d0) go to 62
   55 if (itab(i).ne.0) go to 65
      itab(i)=ifmt(locm+i)
      go to 65
   62 value(locv+i)=defval(locm+i)
   65 continue
      loc=nodplc(loc)
      go to 50
c
c     limit model values
c
   70 go to (80,85,90,95), id
c...  diodes
   80 loc=locate(21)
   82 if (loc.eq.0) go to 130
      locv=nodplc(loc+1)
      value(locv+7)=dmin1(value(locv+7),0.9d0)
      value(locv+8)=dmax1(value(locv+8),0.1d0)
      value(locv+11)=dmax1(value(locv+11),0.1d0)
      value(locv+12)=dmin1(value(locv+12),0.95d0)
      loc=nodplc(loc)
      go to 82
c...  bipolar transistors
   85 loc=locate(22)
   87 if (loc.eq.0) go to 130
      locv=nodplc(loc+1)
      value(locv+23)=dmin1(value(locv+23),0.9d0)
      if (value(locv+24).eq.0.0d0) value(locv+28)=0.0d0
      value(locv+31)=dmin1(value(locv+31),0.9d0)
      value(locv+32)=dmin1(value(locv+32),1.0d0)
      value(locv+40)=dmin1(value(locv+40),0.9d0)
      value(locv+42)=dmax1(value(locv+42),0.1d0)
      value(locv+45)=dmax1(value(locv+45),0.1d0)
      value(locv+46)=dmin1(value(locv+46),0.9999d0)
      loc=nodplc(loc)
      if (value(locv+18).eq.0.0d0) value(locv+18)=value(locv+16)
      if (value(locv+16).ge.value(locv+18)) go to 87
      write(iofile,89) value(locv)
   89 format('0warning:  minimum base resistance (rbm) is less than '
     1       ,'total (rb) for model ',a8,/10x,' rbm set equal to rb',/)
      value(locv+18)=value(locv+16)
      go to 87
c...  jfets
   90 loc=locate(23)
   92 if (loc.eq.0) go to 130
      locv=nodplc(loc+1)
      value(locv+11)=dmax1(value(locv+11),0.1d0)
      value(locv+12)=dmin1(value(locv+12),0.95d0)
      loc=nodplc(loc)
      go to 92
c...  mosfets
   95 loc=locate(24)
   97 if (loc.eq.0) go to 130
      locv=nodplc(loc+1)
c
  100 value(locv+37)=dmax1(value(locv+37),0.1d0)
      value(locv+38)=dmin1(value(locv+38),0.95d0)
      if (value(locv+23).le.0.0d0) go to 120
      cj=dsqrt(epssil*charge*value(locv+23)*1.0d6/
     1   (2.0d0*value(locv+12)))
      if (value(locv+9).le.0.0d0) go to 105
      itab(9)=2
  105 if (value(locv+10).le.0.0d0) go to 110
      itab(10)=2
      go to 115
  110 if (value(locv+17).le.0.0d0) value(locv+17)=cj
      itab(17)=2
  115 if ((value(locv+7).le.0.0d0).and.
     1    (value(locv+8).le.0.0d0)) go to 120
      itab(7)=2
      itab(8)=2
  120 if (value(locv+6).ge.0.2d0) write (iofile,121) value(locv)
  121 format ('0warning:  the value of lambda for mosfet model ',a8,/,
     1   ' is unusually large and might cause nonconvergence',/)
      if (lev.ne.2) value(locv+35)=1.0d0
      if (lev.ne.3) go to 125
      itab(40)=1
      itab(41)=1
      itab(42)=1
      itab(43)=1
  125 loc=nodplc(loc)
      go to 97
c
c     print model parameters
c
  130 if (iprntm.eq.0) go to 360
      locs=locate(id+20)
  140 kntr=0
      loc=locs
      go to (150,160,170,180), id
  150 call title(0,lwidth,1,titled)
      go to 200
  160 call title(0,lwidth,1,titleb)
      go to 200
  170 call title(0,lwidth,1,titlej)
      go to 200
  180 call title(0,lwidth,1,titlem)
  200 if (loc.eq.0) go to 210
      if (kntr.lt.kntlim) go to 220
  210 locn=loc
      go to 240
  220 kntr=kntr+1
      locv=nodplc(loc+1)
      atable(kntr)=value(locv)
  230 loc=nodplc(loc)
      go to 200
  240 write (iofile,241) (atable(k),k=1,kntr)
  241 format(//11x,12(2x,a8))
      if (id.eq.1) go to 300
      kntr=0
      loc=locs
  250 if (loc.eq.0) go to 260
      if (kntr.ge.kntlim) go to 260
      kntr=kntr+1
      atable(kntr)=antype(id)
      if (nodplc(loc+2).eq.-1) atable(kntr)=aptype(id)
      loc=nodplc(loc)
      go to 250
  260 write (iofile,261) (atable(k),k=1,kntr)
  261 format('0type',4x,12(4x,a6))
  300 do 340 i=1,nopar
      if (itab(i).eq.0) go to 340
      kntr=0
      iccflg=0
      loc=locs
  310 if (loc.eq.0) go to 320
      if (kntr.ge.kntlim) go to 320
      locv=nodplc(loc+1)
      kntr=kntr+1
      if (iccflg.ne.0) go to 313
      if (id.ne.2) go to 315
      if ((i.ne.6).and.(i.ne.12)) go to 315
      if (value(locv+i).le.1.0d0) go to 315
      iccflg=i/6
  313 btable(kntr)=value(locv+i)
      value(locv+i)=value(locv+i)*value(locv+1)
  315 atable(kntr)=value(locv+i)
      loc=nodplc(loc)
      go to 310
  320 if (itab(i).eq.2) go to 330
      write (iofile,321) ampar(locm+i),(atable(k),k=1,kntr)
  321 format(1h0,a8,12f10.3)
      go to 340
  330 write (iofile,331) ampar(locm+i),(atable(k),k=1,kntr)
  331 format(1h0,a8,1p12d10.2)
      if (iccflg.eq.0) go to 340
      write (iofile,321) cpar(iccflg),(btable(k),k=1,kntr)
  340 continue
      if (locn.eq.0) go to 390
      locs=locn
      go to 140
c
c  special  treatment for c2 & c4 in the bjt model
c  when no model parameter print
c
  360 if (id.ne.2) go to 390
      loc=locate(id+20)
  370 if (loc.eq.0) go to 390
      locv=nodplc(loc+1)
      if (value(locv+6).ge.1.0d0)
     1   value(locv+6)=value(locv+6)*value(locv+1)
      if (value(locv+12).ge.1.0d0)
     1   value(locv+12)=value(locv+12)*value(locv+1)
      loc=nodplc(loc)
      go to 370
  390 continue
c
c  process model parameters
c
c  diodes
c
  400 loc=locate(21)
  410 if (loc.eq.0) go to 420
      locv=nodplc(loc+1)
      if (value(locv+2).ne.0.0d0) value(locv+2)=1.0d0/value(locv+2)
      pb=value(locv+6)
      xm=value(locv+7)
      fc=value(locv+12)
      value(locv+12)=fc*pb
      xfc=dlog(1.0d0-fc)
      value(locv+15)=pb*(1.0d0-dexp((1.0d0-xm)*xfc))/(1.0d0-xm)
      value(locv+16)=dexp((1.0d0+xm)*xfc)
      value(locv+17)=1.0d0-fc*(1.0d0+xm)
      csat=value(locv+1)
      vte=value(locv+3)*vt
      value(locv+18)=vte*dlog(vte/(root2*csat))
      bv=value(locv+13)
      if (bv.eq.0) go to 418
      cbv=value(locv+14)
      if (cbv.ge.csat*bv/vt) go to 412
      cbv=csat*bv/vt
      write (iofile,411) value(locv),cbv
  411 format('0warning:  in diode model ',a8,' ibv increased to ',1pe10.
     1   3,11x,'to resolve incompatibility with specified is'/)
      xbv=bv
      go to 416
  412 tol=reltol*cbv
      xbv=bv-vt*dlog(1.0d0+cbv/csat)
      iter=0
  413 xbv=bv-vt*dlog(cbv/csat+1.0d0-xbv/vt)
      xcbv=csat*(dexp((bv-xbv)/vt)-1.0d0+xbv/vt)
      if (dabs(xcbv-cbv).le.tol) go to 416
      iter=iter+1
      if (iter.lt.25) go to 413
      write (iofile,415) xbv,xcbv
  415 format('0warning:  unable to match forward and reverse diode regio
     1ns',/,11x,'bv = ',1pd10.3,' and ibv = ',d10.3,/)
  416 value(locv+13)=xbv
  418 loc=nodplc(loc)
      go to 410
c
c  bipolar transistor models
c
  420 loc=locate(22)
  430 if (loc.eq.0) go to 440
      locv=nodplc(loc+1)
      if (value(locv+4).ne.0.0d0) value(locv+4)=1.0d0/value(locv+4)
      if (value(locv+5).ne.0.0d0) value(locv+5)=1.0d0/value(locv+5)
      if (value(locv+10).ne.0.0d0) value(locv+10)=1.0d0/value(locv+10)
      if (value(locv+11).ne.0.0d0) value(locv+11)=1.0d0/value(locv+11)
      if (value(locv+19).ne.0.0d0) value(locv+19)=1.0d0/value(locv+19)
      if (value(locv+20).ne.0.0d0) value(locv+20)=1.0d0/value(locv+20)
      if (value(locv+26).ne.0.0d0) value(locv+26)=1.0d0/value(locv+26)
     1   /1.44d0
      value(locv+28)=value(locv+28)/rad*value(locv+24)
      if (value(locv+35).ne.0.0d0) value(locv+35)=1.0d0/value(locv+35)
     1   /1.44d0
      pe=value(locv+22)
      xme=value(locv+23)
      pc=value(locv+30)
      xmc=value(locv+31)
      fc=value(locv+46)
      value(locv+46)=fc*pe
      xfc=dlog(1.0d0-fc)
      value(locv+47)=pe*(1.0d0-dexp((1.0d0-xme)*xfc))/(1.0d0-xme)
      value(locv+48)=dexp((1.0d0+xme)*xfc)
      value(locv+49)=1.0d0-fc*(1.0d0+xme)
      value(locv+50)=fc*pc
      value(locv+51)=pc*(1.0d0-dexp((1.0d0-xmc)*xfc))/(1.0d0-xmc)
      value(locv+52)=dexp((1.0d0+xmc)*xfc)
      value(locv+53)=1.0d0-fc*(1.0d0+xmc)
      csat=value(locv+1)
      value(locv+54)=vt*dlog(vt/(root2*csat))
      loc=nodplc(loc)
      go to 430
c
c  jfet models
c
  440 loc=locate(23)
  450 if (loc.eq.0) go to 460
      locv=nodplc(loc+1)
      if (value(locv+4).ne.0.0d0) value(locv+4)=1.0d0/value(locv+4)
      if (value(locv+5).ne.0.0d0) value(locv+5)=1.0d0/value(locv+5)
      pb=value(locv+8)
      xm=0.5d0
      fc=value(locv+12)
      value(locv+12)=fc*pb
      xfc=dlog(1.0d0-fc)
      value(locv+13)=pb*(1.0d0-dexp((1.0d0-xm)*xfc))/(1.0d0-xm)
      value(locv+14)=dexp((1.0d0+xm)*xfc)
      value(locv+15)=1.0d0-fc*(1.0d0+xm)
      csat=value(locv+9)
      value(locv+16)=vt*dlog(vt/(root2*csat))
      loc=nodplc(loc)
      go to 450
c
c  mosfet models
c
  460 loc=locate(24)
  470 if (loc.eq.0) go to 600
      locv=nodplc(loc+1)
      type=nodplc(loc+2)
      if (value(locv+7).ne.0.0d0) value(locv+7)=1.0d0/value(locv+7)
      if (value(locv+8).ne.0.0d0) value(locv+8)=1.0d0/value(locv+8)
      if (value(locv+16).ne.0.0d0) value(locv+16)=1.0d0/value(locv+16)
      value(locv+23)=value(locv+23)*1.0d6
      value(locv+24)=value(locv+24)*1.0d4
      value(locv+25)=value(locv+25)*1.0d4
      if (value(locv+22).ne.0.0d0) value(locv+22)=epsox/value(locv+22)
      value(locv+29)=value(locv+29)*1.0d-4
      if (lev.eq.3) go to 472
      value(locv+30)=value(locv+30)*1.0d2
      go to 473
c
c   move mos3 parameters : theta from locations locv+40 to locv+30
c                          eta                       41         31
c                          kappa                     42         32
c   and replace locv+6 by (xd)**2
c
  472 value(locv+39)=value(locv+39)
     1   *0.25d0*twopi*epssil/value(locv+22)
      value(locv+30)=value(locv+40)
      value(locv+31)=value(locv+41)*8.15d-22/value(locv+22)
      value(locv+32)=value(locv+42)
      if (value(locv+23).gt.0.0d0)
     1    value(locv+6)=(epssil+epssil)/(charge*value(locv+23))
c
c   noise parameters
c
  473 pb=value(locv+12)
      xm=0.5d0
      fc=value(locv+38)
      value(locv+38)=fc*pb
      xfc=dlog(1.0d0-fc)
      value(locv+40)=pb*(1.0d0-dexp((1.0d0-xm)*xfc))/(1.0d0-xm)
      value(locv+41)=dexp((1.0d0+xm)*xfc)
      value(locv+42)=1.0d0-fc*(1.0d0+xm)
      value(locv+43)=-1.0d0
      value(locv+44)=value(locv+2)-
     1   type*value(locv+4)*dsqrt(value(locv+5))
  475 if (value(locv+22).ne.0.0d0.and.lev.ne.3)
     1   value(locv+30)=value(locv+30)*epssil/value(locv+22)
      loc=nodplc(loc)
      go to 470
c
c  reserve additional nodes
c
c  diodes
c
  600 loc=locate(11)
  610 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 700
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      if (value(locm+2).eq.0.0d0) go to 620
      numnod=numnod+1
      nodplc(loc+4)=numnod
      go to 630
  620 nodplc(loc+4)=nodplc(loc+2)
  630 loc=nodplc(loc)
      go to 610
c
c  transistors
c
  700 loc=locate(12)
  710 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 800
      nodplc(loc+30)=nodplc(loc+5)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      if (value(locm+16).eq.0.0d0) go to 720
      numnod=numnod+1
      nodplc(loc+6)=numnod
      go to 730
  720 nodplc(loc+6)=nodplc(loc+3)
  730 if (value(locm+20).eq.0.0d0) go to 740
      numnod=numnod+1
      nodplc(loc+5)=numnod
      go to 750
  740 nodplc(loc+5)=nodplc(loc+2)
  750 if (value(locm+19).eq.0.0d0) go to 760
      numnod=numnod+1
      nodplc(loc+7)=numnod
      go to 770
  760 nodplc(loc+7)=nodplc(loc+4)
  770 loc=nodplc(loc)
      go to 710
c
c  jfets
c
  800 loc=locate(13)
  810 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 900
      locm=nodplc(loc+7)
      locm=nodplc(locm+1)
      if (value(locm+4).eq.0.0d0) go to 820
      numnod=numnod+1
      nodplc(loc+5)=numnod
      go to 830
  820 nodplc(loc+5)=nodplc(loc+2)
  830 if (value(locm+5).eq.0.0d0) go to 840
      numnod=numnod+1
      nodplc(loc+6)=numnod
      go to 850
  840 nodplc(loc+6)=nodplc(loc+4)
  850 loc=nodplc(loc)
      go to 810
c
c  mosfets
c
  900 loc=locate(14)
  910 if (loc.eq.0) go to 1000
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      locv=nodplc(loc+1)
      xleff=value(locv+1)-2.0d0*value(locm+28)
      if (xleff.gt.0.0d0) go to 915
      write(iofile,911) value(locv),value(locm)
  911 format('0*error*:  effective channel length of ',a8,' less than ',
     1   'zero.',/' check value of ld for model ',a8)
      if (nodplc(loc+33).ne.0) go to 960
  915 if ((value(locm+7).eq.0.0d0).and.
     1    (value(locm+16).eq.0.0d0)) go to 920
      numnod=numnod+1
      nodplc(loc+6)=numnod
      go to 930
  920 nodplc(loc+6)=nodplc(loc+2)
  930 if ((value(locm+8).eq.0.0d0).and.
     1    (value(locm+16).eq.0.0d0)) go to 940
      numnod=numnod+1
      nodplc(loc+7)=numnod
      go to 950
  940 nodplc(loc+7)=nodplc(loc+4)
  950 ad=value(locv+3)
      as=value(locv+4)
      if ((ad.le.0.0d0).or.(as.le.0.0d0)
     1   .and.value(locm+11).le.0.0d0)
     2   value(locm+11)=1.0d-14
  960 loc=nodplc(loc)
      go to 910
c
c  transmission lines
c
 1000 loc=locate(17)
 1010 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 2000
      numnod=numnod+1
      nodplc(loc+6)=numnod
      numnod=numnod+1
      nodplc(loc+7)=numnod
      loc=nodplc(loc)
      go to 1010
c
c  finished
c
 2000 return
      end
