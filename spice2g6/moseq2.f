      subroutine moseq2(vds,vbs,vgs,gm,gds,gmbs,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      implicit double precision (a-h,o-z)
c
c     this routine evaluates the drain current, its derivatives and
c     the charges associated with the gate, channel and bulk
c     for mosfets
c
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
c
      dimension a4(4),b4(4),x4(8),poly4(8),sig1(4),sig2(4)
      data sig1 / 1.0d0, -1.0d0, 1.0d0, -1.0d0/,
     1     sig2 / 1.0d0,  1.0d0,-1.0d0, -1.0d0/
c
c     icharg=1 causes charges to be computed
c     icharg=0 bypasses the computation of charges
c
      icharg=1
      if (mode.ne.1.and.xqco.le.0.5d0) go to 100
      icharg=0
      if (xqco.gt.0.5d0) go to 100
      if (modedc.eq.2.and.nosolv.ne.0) icharg=1
      if (initf.eq.4) icharg=1
c
c  compute some useful quantities
c
  100 if (vbs.gt.0.0d0) go to 110
      sarg=dsqrt(phi-vbs)
      tsarg=sarg+sarg
      dsrgdb=-0.5d0/sarg
      d2sdb2=+0.5d0*dsrgdb/(phi-vbs)
      go to 120
  110 sphi=dsqrt(phi)
      sphi3=phi*sphi
      sarg=sphi/(1.0d0+0.5d0*vbs/phi)
      tsarg=sarg+sarg
      dsrgdb=-0.5d0*sarg*sarg/sphi3
      d2sdb2=-dsrgdb*sarg/sphi3
  120 if ((vds-vbs).lt.0.0d0) go to 130
      barg=dsqrt(phi+vds-vbs)
      dbrgdb=-0.5d0/barg
      d2bdb2=+0.5d0*dbrgdb/(phi+vds-vbs)
      go to 200
  130 barg=sphi/(1.0d0+0.5d0*(vbs-vds)/phi)
      dbrgdb=-0.5d0*barg*barg/sphi3
      d2bdb2=-dbrgdb*barg/sphi3
c
c  calculate threshold voltage (von)
c     narrow-channel effect
c
  200 factor=0.125d0*fnarrw*twopi*epssil/cox*xl
      eta=1.0d0+factor
      vbin=vbi+factor*(phi-vbs)
      if (gamma.le.0.0d0) go to 215
      if (xnsub.le.0.0d0) go to 215
      xwd=xd*barg
      xws=xd*sarg
c
c     short-channel effect with vds .ne. 0.0d0
c
      argss=0.0d0
      argsd=0.0d0
      dbargs=0.0d0
      dbargd=0.0d0
      dgdvds=0.0d0
      dgddb2=0.0d0
      if (xj.le.0.0d0) go to 205
      argxs=1.0d0+2.0d0*xws/xj
      args=dsqrt(argxs)
      argss=0.5d0*xj/xl*(args-1.0d0)
      argxd=1.0d0+2.0d0*xwd/xj
      argd=dsqrt(argxd)
      argsd=0.5d0*xj/xl*(argd-1.0d0)
  205 gamasd=gamma*(1.0d0-argss-argsd)
      gamass=gamma*(1.0d0-2.0d0*argss)
      dbxwd=xd*dbrgdb
      dbxws=xd*dsrgdb
      if (xj.le.0.0d0) go to 210
      dbargs=0.5d0/xl*dbxws/args
      dbargd=0.5d0/xl*dbxwd/argd
      dasdb2=-xd*( d2sdb2+dsrgdb*dsrgdb*xd/(xj*argxs) )/(xl*args)
      daddb2=-xd*( d2bdb2+dbrgdb*dbrgdb*xd/(xj*argxd) )/(xl*argd)
      dgddb2=-0.5d0*gamma*(dasdb2+daddb2)
  210 dgddvb=-gamma*(dbargs+dbargd)
      dgsdvb=-2.0d0*gamma*dbargs
      if (xj.le.0.0d0) go to 220
      ddxwd=-dbxwd
      dgdvds=-gamma*0.5d0/xl*ddxwd/argd
      go to 220
  215 gamasd=gamma
      gamass=gamma
      gammad=gamma
      dgddvb=0.0d0
      dgsdvb=0.0d0
      dgdvds=0.0d0
      dgddb2=0.0d0
  220 von=vbin+gamasd*sarg
c     write(iofile,221) von,vbin,vbi,gamasd,argss,argsd,xj
  221 format ('0msg1:'/1p7d10.2)
      vth=von
      vdsat=0.0d0
  225 if (xnfs.eq.0.0d0.or.cox.eq.0.0d0) go to 230
      cfs=charge*xnfs
      cdonco=-(gamasd*dsrgdb+dgddvb*sarg)+factor
      xn=1.0d0+cfs/cox*xw*xl+cdonco
      von=von+vt*xn
c     write (iofile,226) von,cdonco,xn,cfs,xd
  226 format(' msg2:'/1p6d10.2)
      argg=1.0d0/(vt*xn)
      vgst=vgs-von
      go to 300
  230 vgst=vgs-von
      if (vgs.gt.von) go to 300
c
c  cutoff region
c
      gds=0.0d0
      go to 1050
c
c  compute some more useful quantities
c
  300 sarg3=sarg*sarg*sarg
      sbiarg=dsqrt(phib)
      gammad=gamasd
      dgdvbs=dgddvb
      body=barg*barg*barg-sarg3
      gdbdv=2.0d0*gammad*(barg*barg*dbrgdb-sarg*sarg*dsrgdb)
      dodvbs=-factor+dgdvbs*sarg+gammad*dsrgdb
      if (xnfs.eq.0.0d0) go to 400
      if (cox.eq.0.0d0) go to 410
      dxndvb=2.0d0*dgdvbs*dsrgdb+gammad*d2sdb2+dgddb2*sarg
      dodvbs=dodvbs+vt*dxndvb
      dxndvd=dgdvds*dsrgdb
      dodvds=dgdvds*sarg+vt*dxndvd
c
c  evaluate effective mobility and its derivatives
c
  400 if (cox.le.0.0d0) go to 410
      udenom=vgst
      if (udenom.le.vbp) go to 410
      ufact=dexp(uexp*dlog(vbp/udenom))
      ueff=uo*ufact
      dudvgs=-ufact*uexp/udenom
      dudvds=0.0d0
      dudvbs=uexp*ufact*dodvbs/vgst
      go to 500
  410 ufact=1.0d0
      ueff=uo
      dudvgs=0.0d0
      dudvds=0.0d0
      dudvbs=0.0d0
c
c     evaluate saturation voltage and its derivatives according to
c     grove-frohman equation
c
  500 vgsx=vgs
      gammad=gamasd/eta
      dgdvbs=dgddvb
      if (xnfs.ne.0.0d0.and.cox.ne.0.0d0)
     1   vgsx=dmax1(vgs,von)
  505 if (gammad.le.0.0d0) go to 535
      gammd2=gammad*gammad
      argv=(vgsx-vbin)/eta+phi-vbs
      if (argv.le.0.0d0) go to 540
      arg=dsqrt(1.0d0+4.0d0*argv/gammd2)
      vdsat=(vgsx-vbin)/eta+gammd2*(1.0d0-arg)/2.0d0
      vdsat=dmax1(vdsat,0.0d0)
  510 if (icharg.eq.0) go to 530
      arg1=gammd2/(eta*eta)
      arg2=vds-0.5d0*arg1
      argsq=(arg2+0.5d0*arg1+phi-vbs)*arg1
      if (argsq.ge.0.0d0) go to 515
      vpof=vth
      go to 520
  515 vpof=vbin+eta*(arg2+0.5d0*arg1+dsqrt(argsq))
  520 argv1=(vpof-vbin)/eta+phi-vbs
      if (argv1.gt.0.0d0) go to 525
      vdsat1=0.0d0
      go to 530
  525 arg1=dsqrt(1.0d0+4.0d0*argv1/gammd2)
      vdsat1=(vpof-vbin)/eta+gammd2*(1.0d0-arg1)/2.0d0
      vdsat1=dmax1(vdsat1,0.0d0)
  530 dsdvgs=(1.0d0-1.0d0/arg)/eta
      dsdvbs=(gammad*(1.0d0-arg)+2.0d0*argv/(gammad*arg))/eta*dgdvbs+
     1       1.0d0/arg+factor*dsdvgs
      go to 545
  535 vdsat=dmax1((vgsx-vbin)/eta,0.0d0)
      vpof=dmax1((eta*vds+vbin),0.0d0)
      vdsat1=dmax1((vpof-vbin)/eta,0.0d0)
      dsdvgs=1.0d0
      dsdvbs=0.0d0
      go to 545
  540 vdsat=0.0d0
      vpof=vth
      vdsat1=0.0d0
      dsdvgs=0.0d0
      dsdvbs=0.0d0
c
c     store vdsat as above in vpof (pinch-off)
c
  545 if (vmax.le.0.0d0) go to 600
c
c     evaluate saturation voltage and its derivatives according to
c     baum's theory of scattering velocity saturation
c
      gammd2=gammad*gammad
      v1=(vgsx-vbin)/eta+phi-vbs
      v2=phi-vbs
      xv=vmax*xl/ueff
      a1=gammad/0.75d0
      b1=-2.0d0*(v1+xv)
      c1=-2.0d0*gammad*xv
      d1=2.0d0*v1*(v2+xv)-v2*v2-4.0d0/3.0d0*gammad*sarg3
      a=-b1
      b=a1*c1-4.0d0*d1
      c=-d1*(a1*a1-4.0d0*b1)-c1*c1
      r=-a*a/3.0d0+b
      s=2.0d0*a*a*a/27.0d0-a*b/3.0d0+c
      r3=r*r*r
      s2=s*s
      p=s2/4.0d0+r3/27.0d0
      p0=dabs(p)
      p2=dsqrt(p0)
      if (p.ge.0.0d0) go to 550
      ro=dsqrt(s2/4.0d0+p0)
      ro=dlog(ro)/3.0d0
      ro=dexp(ro)
      fi=datan(-2.0d0*p2/s)
      y3=2.0d0*ro*dcos(fi/3.0d0)-a/3.0d0
      go to 560
  550 p3=dexp(dlog(dabs(-s/2.0d0+p2))/3.0d0)
      p4=dexp(dlog(dabs(-s/2.0d0-p2))/3.0d0)
      y3=p3+p4-a/3.0d0
  560 iknt=0
      a3=dsqrt(a1*a1/4.0d0-b1+y3)
      b3=dsqrt(y3*y3/4.0d0-d1)
      do 570 i=1,4
      a4(i)=a1/2.0d0+sig1(i)*a3
      b4(i)=y3/2.0d0+sig2(i)*b3
      delta4=a4(i)*a4(i)/4.0d0-b4(i)
      if (delta4.lt.0.0d0) go to 570
      iknt=iknt+1
      x4(iknt)=-a4(i)/2.0d0+dsqrt(delta4)
      iknt=iknt+1
      x4(iknt)=-a4(i)/2.0d0-dsqrt(delta4)
  570 continue
      jknt=0
      do 580 j=1,iknt
      if (x4(j).le.0.0d0) go to 580
      poly4(j)=x4(j)*x4(j)*x4(j)*x4(j)+a1*x4(j)*x4(j)*x4(j)
      poly4(j)=poly4(j)+b1*x4(j)*x4(j)+c1*x4(j)+d1
      if (dabs(poly4(j)).gt.1.0d-6) go to 580
      jknt=jknt+1
      if (jknt.gt.1) go to 575
      xvalid=x4(j)
  575 if (x4(j).gt.xvalid) go to 580
      xvalid=x4(j)
  580 continue
      if (jknt.gt.0) go to 590
      ivmflg=ivmflg+1
      go to 600
  590 vdsat=xvalid*xvalid+vbs-phi
c
c  evaluate effective channel length and its derivatives
c
  600 if (vds.eq.0.0d0) go to 610
      gammad=gamasd
      if ((vbs-vdsat).gt.0.0d0) go to 601
      bsarg=dsqrt(vdsat-vbs+phi)
      dbsrdb=-0.5d0/bsarg
      go to 602
  601 bsarg=sphi/(1.0d0+0.5d0*(vbs-vdsat)/phi)
      dbsrdb=-0.5d0*bsarg*bsarg/sphi3
  602 bodys=bsarg*bsarg*bsarg-sarg3
      gdbdvs=2.0d0*gammad*(bsarg*bsarg*dbsrdb-sarg*sarg*dsrgdb)
      if (vmax.gt.0.0d0) go to 603
      if (xnsub.eq.0.0d0) go to 610
      if (xlamda.gt.0.0d0) go to 610
      argv=(vds-vdsat)/4.0d0
      sargv=dsqrt(1.0d0+argv*argv)
      arg=dsqrt(argv+sargv)
      xlfact=xd/(xl*vds)
      xlamda=xlfact*arg
      dldsat=vds*xlfact*arg/(8.0d0*sargv)
      go to 605
  603 argv=(vgsx-vbin)/eta-vdsat
      xdv=xd/dsqrt(xneff)
      xlv=vmax*xdv/(2.0d0*ueff)
      vqchan=argv-gammad*bsarg
      dqdsat=-1.0d0+gammad*dbsrdb
      vl=vmax*xl
      dfunds=vl*dqdsat-ueff*vqchan
      dfundg=(vl-ueff*vdsat)/eta
      dfundb=-vl*(1.0d0+dqdsat-factor/eta)+
     1        ueff*(gdbdvs-dgdvbs*bodys/1.5d0)/eta
      dsdvgs=-dfundg/dfunds
      dsdvbs=-dfundb/dfunds
      if (xnsub.eq.0.0d0) go to 610
      if (xlamda.gt.0.0d0) go to 610
      argv=dmax1(vds-vdsat,0.0d0)
      xls=dsqrt(xlv*xlv+argv)
      dldsat=xdv/(2.0d0*xls)
      xlfact=xdv/(xl*vds)
      xlamda=xlfact*(xls-xlv)
      dldsat=dldsat/xl
  605 dldvgs=dldsat*dsdvgs
      dldvds=-xlamda+dldsat
      dldvbs=dldsat*dsdvbs
      go to 620
  610 dldvgs=0.0d0
      dldvds=0.0d0
      dldvbs=0.0d0
c
c     limit channel shortening at punch-through
c
  620 xwb=xd*sbiarg
      xld=xl-xwb
      clfact=1.0d0-xlamda*vds
      dldvds=-xlamda-dldvds
      xleff=xl*clfact
      deltal=xlamda*vds*xl
      if (xnsub.eq.0.0d0) xwb=0.25d-6
      if (xleff.ge.xwb) go to 700
      xleff=xwb/(1.0d0+(deltal-xld)/xwb)
      clfact=xleff/xl
      dfact=xleff*xleff/(xwb*xwb)
      dldvgs=dfact*dldvgs
      dldvds=dfact*dldvds
      dldvbs=dfact*dldvbs
c
c  evaluate effective beta (effective kp)
c
  700 beta1=beta*ufact/clfact
c
c  test for mode of operation and branch appropriately
c
      gammad=gamasd
      dgdvbs=dgddvb
      if (vds.gt.1.0d-8) go to 730
      if (vgs.gt.von) go to 720
      if ((xnfs.ne.0.0d0).and.(cox.ne.0.0d0)) go to 710
      gds=0.0d0
      go to 1050
c
  710 gds=beta1*(von-vbin-gammad*sarg)*dexp(argg*(vgs-von))
      go to 1050
c
c
  720 gds=beta1*(vgs-vbin-gammad*sarg)
      go to 1050
c
  730 if (vgs.gt.von) go to 900
c
c  subthreshold region
c
      if (vdsat.gt.0.0d0) go to 830
      gds=0.0d0
      if (vgs.gt.vth) go to 1020
      go to 1050
  830 vdson=dmin1(vdsat,vds)
      if (vds.le.vdsat) go to 850
      barg=bsarg
      dbrgdb=dbsrdb
      body=bodys
      gdbdv=gdbdvs
  850 cdson=beta1*((von-vbin-eta*vdson*0.5d0)*vdson-gammad*body/1.5d0)
      didvds=beta1*(von-vbin-eta*vdson-gammad*barg)
      gdson=-cdson*dldvds/clfact-beta1*dgdvds*body/1.5d0
      if (vds.lt.vdsat) gdson=gdson+didvds
      gbson=-cdson*dldvbs/clfact
     1      +beta1*(dodvbs*vdson+factor*vdson-dgdvbs*body/1.5d0-gdbdv)
      if (vds.gt.vdsat) gbson=gbson+didvds*dsdvbs
      expg=dexp(argg*(vgs-von))
      cdrain=cdson*expg
      gmw=cdrain*argg
      gm=gmw
      if (vds.gt.vdsat) gm=gmw+didvds*dsdvgs*expg
      gds=gdson*expg-gm*dodvds-gmw*(vgs-von)*dxndvd/xn
      gmbs=gbson*expg-gm*dodvbs-gmw*(vgs-von)*dxndvb/xn
      go to 1020
c
c
  900 if (vds.gt.vdsat) go to 1000
c
c  linear region
c
      cdrain=beta1*((vgs-vbin-eta*vds/2.0d0)*vds-gammad*body/1.5d0)
      arg=cdrain*(dudvgs/ufact-dldvgs/clfact)
      gm=arg+beta1*vds
      arg=cdrain*(dudvds/ufact-dldvds/clfact)
      gds=arg+beta1*(vgs-vbin-eta*vds-
     1   gammad*barg-dgdvds*body/1.5d0)
      arg=cdrain*(dudvbs/ufact-dldvbs/clfact)
      gmbs=arg-beta1*(gdbdv+dgdvbs*body/1.5d0-factor*vds)
      go to 1020
c
c  saturation region
c
 1000 cdrain=beta1*((vgs-vbin-eta*vdsat/2.0d0)*vdsat-gammad*bodys/1.5d0)
      arg=cdrain*(dudvgs/ufact-dldvgs/clfact)
      gm=arg+beta1*vdsat+
     1   beta1*(vgs-vbin-eta*vdsat-gammad*bsarg)*dsdvgs
      gds=-cdrain*dldvds/clfact-beta1*dgdvds*bodys/1.5d0
      arg=cdrain*(dudvbs/ufact-dldvbs/clfact)
      gmbs=arg-beta1*(gdbdvs+dgdvbs*bodys/1.5d0-factor*vdsat)+
     1     beta1*(vgs-vbin-eta*vdsat-gammad*bsarg)*dsdvbs
c
c     compute charges for "on" region
c
 1020 if (icharg.eq.0) go to 1500
      if (vgs.le.vth) go to 1070
      call mqspof(vds,vbs,vgs,vpof,vdsat1,vth,vbin,gamasd,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      go to 2000
c
c  finish special cases
c
 1050 cdrain=0.0d0
      gm=0.0d0
      gmbs=0.0d0
 1070 xqc=xqco
      if (icharg.eq.0) go to 1500
      call mosq2(vds,vbs,vgs,vdsat,vth,vbin,gamasd,cox,phi,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      qspof=0.0d0
      go to 2000
c
c  finished
c
 1500 qg=0.0d0
      qb=0.0d0
      qc=0.0d0
      qspof=0.0d0
 2000 return
      end
