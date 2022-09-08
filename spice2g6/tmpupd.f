      subroutine tmpupd
      implicit double precision (a-h,o-z)
c
c     this routine updates the temperature-dependent parameters in the
c device models.  it also updates the values of temperature-dependent
c resistors.  the updated values are printed.
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
      dimension tmptit(4)
      data tmptit / 8htemperat, 8hure-adju, 8hsted val, 8hues      /
c
c
      reftmp=27.0d0+ctok
      temp1=value(itemps+itemno-1)+ctok
      temp2=value(itemps+itemno)+ctok
      xkt=boltz*temp2
      oldvt=vt
      vt=xkt/charge
      oldeg=egfet
      egfet=1.16d0-(7.02d-4*temp2*temp2)/(temp2+1108.0d0)
      arg=-egfet/(xkt+xkt)+1.1150877d0/(boltz*(reftmp+reftmp))
      ratio=temp2/temp1
      ratlog=dlog(ratio)
      ratio1=ratio-1.0d0
      dtemp=temp2-reftmp
      delt=value(itemps+itemno)-value(itemps+1)
      deltsq=delt*delt
      fact2=temp2/reftmp
      xni=1.45d16*fact2*dsqrt(fact2)*dexp(charge*arg)
      pbfact=-2*vt*(1.5d0*dlog(fact2)+charge*arg)
      xkt1=boltz*temp1
      vt1=xkt1/charge
      egfet1=1.16d0-(7.02d-4*temp1*temp1)/(temp1+1108.0d0)
      arg1=-egfet1/(xkt1+xkt1)+1.1150877d0/(boltz*(reftmp+reftmp))
      fact1=temp1/reftmp
      pbfat1=-2*vt1*(1.5d0*dlog(fact1)+charge*arg1)
    5 call title(0,lwidth,1,tmptit)
c
c  resistors
c
      loc=locate(1)
      ititle=0
   10 if (loc.eq.0) go to 100
      locv=nodplc(loc+1)
      tc1=value(locv+3)
      tc2=value(locv+4)
      if (tc1.ne.0.0d0) go to 20
      if (tc2.eq.0.0d0) go to 40
   20 if (ititle.ne.0) go to 30
      write (iofile,21)
   21 format(//'0**** resistors',/,'0name',8x,'value',//)
      ititle=1
   30 rnew=value(locv+2)*(1.0d0+tc1*delt+tc2*deltsq)
      value(locv+1)=1.0d0/rnew
      write (iofile,31) value(locv),rnew
   31 format(1x,a8,1p6d11.3)
   40 loc=nodplc(loc)
      go to 10
c
c  diode model
c
  100 loc=locate(21)
      if (loc.eq.0) go to 200
      write (iofile,101)
  101 format(//'0**** diode model parameters',/,'0name',9x,'is',9x,'vj',
     1   8x,'cjo',//)
  110 if (loc.eq.0) go to 200
      locv=nodplc(loc+1)
c...  is(t2)=is(t1)*dexp(eg/(n*vt)*(t2/t1-1))*(t2/t1)**(xti/n)
      xn=value(locv+3)
      factor=ratio1*value(locv+8)/(xn*vt)+value(locv+9)/xn*ratlog
      factor=dexp(factor)
      value(locv+1)=value(locv+1)*factor
      oldpb=value(locv+6)
      pbo=(value(locv+6)-pbfat1)/fact1
      gmaold=(oldpb-pbo)/pbo
      value(locv+5)=value(locv+5)/(1.0d0+value(locv+7)
     1     *(400.0d-6*(temp1-reftmp)-gmaold))
  120 value(locv+6)=fact2*pbo+pbfact
      gmanew=(value(locv+6)-pbo)/pbo
      value(locv+5)=value(locv+5)
     1   *(1.0d0+value(locv+7)*(400.0d-6*dtemp-gmanew))
      pbrat=value(locv+6)/oldpb
      value(locv+12)=value(locv+12)*pbrat
      value(locv+15)=value(locv+15)*pbrat
      vte=value(locv+3)*vt
      value(locv+18)=vte*dlog(vte/(root2*value(locv+1)))
      write (iofile,31) value(locv),value(locv+1),value(locv+6),
     1                  value(locv+5)
      loc=nodplc(loc)
      go to 110
c
c  bipolar transistor model
c
  200 loc=locate(22)
      if (loc.eq.0) go to 300
      write (iofile,201)
  201 format(//'0**** bjt model parameters',/,'0name',9x,'js',8x,'bf ',
     1   7x,'ise',7x,'br ',7x,'isc',7x,'vje',7x,'cje',7x,'vjc',
     2   7x,'cjc',//)
  210 if (loc.eq.0) go to 300
      locv=nodplc(loc+1)
c...  is(t2)=is(t1)*dexp(eg/vt*(t2/t1-1))*(t2/t1)**xti
      factln=ratio1*value(locv+42)/vt+value(locv+43)*ratlog
      factor=dexp(factln)
      value(locv+1)=value(locv+1)*factor
      tb=value(locv+41)
      bfactr=dexp(tb*ratlog)
      value(locv+2)=value(locv+2)*bfactr
      value(locv+8)=value(locv+8)*bfactr
      value(locv+6)=value(locv+6)*dexp(factln/value(locv+7))/bfactr
      value(locv+12)=value(locv+12)*dexp(factln/value(locv+13))
     1               /bfactr
      oldpb=value(locv+22)
      pbo=(value(locv+22)-pbfat1)/fact1
      gmaold=(oldpb-pbo)/pbo
      value(locv+21)=value(locv+21)/(1.0d0+value(locv+23)
     1     *(400.0d-6*(temp1-reftmp)-gmaold))
  220 value(locv+22)=fact2*pbo+pbfact
      gmanew=(value(locv+22)-pbo)/pbo
      value(locv+21)=value(locv+21)
     1   *(1.0d0+value(locv+23)*(400.0d-6*dtemp-gmanew))
      pbrat=value(locv+22)/oldpb
      value(locv+46)=value(locv+46)*pbrat
      value(locv+47)=value(locv+47)*pbrat
      oldpb=value(locv+30)
      pbo=(value(locv+30)-pbfat1)/fact1
      gmaold=(oldpb-pbo)/pbo
      value(locv+29)=value(locv+29)/(1.0d0+value(locv+31)
     1     *(400.0d-6*(temp1-reftmp)-gmaold))
  230 value(locv+30)=fact2*pbo+pbfact
      gmanew=(value(locv+30)-pbo)/pbo
      value(locv+29)=value(locv+29)
     1   *(1.0d0+value(locv+31)*(400.0d-6*dtemp-gmanew))
      pbrat=value(locv+30)/oldpb
      value(locv+50)=value(locv+50)*pbrat
      value(locv+51)=value(locv+51)*pbrat
      value(locv+54)=vt*dlog(vt/(root2*value(locv+1)))
      write (iofile,211) value(locv),value(locv+1),value(locv+2),
     1   value(locv+6),value(locv+8),value(locv+12),value(locv+22),
     2   value(locv+21),value(locv+30),value(locv+29)
  211 format(1x,a8,1p9d10.3)
      loc=nodplc(loc)
      go to 210
c
c  jfet model
c
  300 loc=locate(23)
      if (loc.eq.0) go to 400
      write (iofile,301)
  301 format(//'0**** jfet model parameters',/,'0name',9x,'is',9x,'pb',
     1   8x,'cgs',8x,'cgd',//)
  310 if (loc.eq.0) go to 400
      locv=nodplc(loc+1)
      value(locv+9)=value(locv+9)*dexp(ratio1*1.11d0/vt)
      oldpb=value(locv+8)
      pbo=(value(locv+8)-pbfat1)/fact1
      gmaold=(oldpb-pbo)/pbo
      oldcjf=1.0d0+0.5d0*(400.0d-6*(temp1-reftmp)-gmaold)
      value(locv+6)=value(locv+6)/oldcjf
      value(locv+7)=value(locv+7)/oldcjf
  320 value(locv+8)=fact2*pbo+pbfact
      gmanew=(value(locv+8)-pbo)/pbo
      cjfact=1.0d0+0.5d0*(400.0d-6*dtemp-gmanew)
      value(locv+6)=value(locv+6)*cjfact
      value(locv+7)=value(locv+7)*cjfact
      pbrat=value(locv+8)/oldpb
      value(locv+12)=value(locv+12)*pbrat
      value(locv+13)=value(locv+13)*pbrat
      value(locv+16)=vt*dlog(vt/(root2*value(locv+9)))
      write (iofile,31) value(locv),value(locv+9),value(locv+8),
     1   value(locv+6),value(locv+7)
      loc=nodplc(loc)
      go to 310
c
c  mosfet model
c
  400 loc=locate(24)
      iprnt=1
  410 if (loc.eq.0) go to 1000
      locv=nodplc(loc+1)
      type=nodplc(loc+2)
      if(iprnt.ne.0) write (iofile,401)
  401 format(//'0**** mosfet model parameters',/,'0name',8x,'vto',8x,
     1   'phi',9x,'pb',7x,'is(js)',7x,'kp',9x,'uo'//)
      iprnt=0
      ratio4=ratio*dsqrt(ratio)
      value(locv+3)=value(locv+3)/ratio4
      value(locv+29)=value(locv+29)/ratio4
      oldphi=value(locv+5)
      phio=(value(locv+5)-pbfat1)/fact1
  415 value(locv+5)=fact2*phio+pbfact
      phi=value(locv+5)
      vfb=value(locv+44)-type*0.5d0*oldphi
      vfb=vfb+0.5d0*(oldeg-egfet)
      value(locv+44)=vfb+type*0.5d0*phi
      value(locv+2)=value(locv+44)+type*value(locv+4)*dsqrt(phi)
      value(locv+11)=value(locv+11)*dexp(-egfet/vt+oldeg/oldvt)
      value(locv+21)=value(locv+21)*dexp(-egfet/vt+oldeg/oldvt)
      oldpb=value(locv+12)
      pbo=(value(locv+12)-pbfat1)/fact1
      gmaold=(oldpb-pbo)/pbo
      coeold=1.0d0+value(locv+18)*(400.0d-6*(temp1-reftmp)-gmaold)
      value(locv+9)=value(locv+9)/coeold
      value(locv+10)=value(locv+10)/coeold
      value(locv+17)=value(locv+17)/coeold
      value(locv+19)=value(locv+19)/(1.0d0+value(locv+20)
     1     *(400.0d-6*(temp1-reftmp)-gmaold))
  420 value(locv+12)=fact2*pbo+pbfact
      gmanew=(value(locv+12)-pbo)/pbo
      coenew=1.0d0+value(locv+18)*(400.0d-6*dtemp-gmanew)
      value(locv+9)=value(locv+9)*coenew
      value(locv+10)=value(locv+10)*coenew
      value(locv+17)=value(locv+17)*coenew
      value(locv+19)=value(locv+19)*
     1   (1.0d0+value(locv+20)*(400.0d-6*dtemp-gmanew))
      pbrat=value(locv+12)/oldpb
      value(locv+37)=value(locv+37)*pbrat
      value(locv+38)=value(locv+38)*pbrat
      csat=dmax1(value(locv+11),value(locv+21))
      write (iofile,31) value(locv),value(locv+2),value(locv+5),
     1   value(locv+12),csat,value(locv+3),value(locv+29)
  430 loc=nodplc(loc)
      go to 410
c
c  finished
c
 1000 return
      end
