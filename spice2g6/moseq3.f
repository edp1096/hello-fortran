      subroutine moseq3(vds,vbs,vgs,gm,gds,gmbs,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      implicit double precision (a-h,o-z)
c
c     this routine evaluates the drain current, its derivatives and
c     the charges associated with the gate, channel and bulk
c     for mosfets based on semi-empirical equations
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
      equivalence (xlamda,alpha),(vbp,theta),(uexp,eta),(utra,xkappa)
      data coeff0/0.0631353d0/,coeff1/0.8013292d0/,coeff2/-0.01110777d0/
c
c     icharg=1 causes charges to be computed
c     icharg=0 bypasses the computation of charges
c
c     icharg=1
c     if (mode.ne.1) go to 10
c     icharg=0
c     if (modedc.eq.2.and.nosolv.ne.0) icharg=1
c     if (initf.eq.4) icharg=1
c
c     reference cdrain equations to source and
c     charge equations to bulk
c
10    continue
      icharg=0
      vgb=vgs-vbs
      vfb=vbi-phi
      vdsat=0.0d0
      qg=0.0d0
      qb=0.0d0
      qc=0.0d0
      cgdb=0.0d0
      cbdb=0.0d0
      onxl=1.0d0/xl
      eta=eta/(xl*xl*xl)
c
c.....square root term
c
      if ( vbs.gt.0.0d0 ) go to 120
      phibs=phi-vbs
      sqphbs=dsqrt(phibs)
      dsqdvb=-0.5d0/sqphbs
      go to 200
120   continue
      sqphis=dsqrt(phi)
      sqphs3=phi*sqphis
      sqphbs=sqphis/(1.0d0+vbs/(phi+phi))
      phibs=sqphbs*sqphbs
      dsqdvb=-phibs/(sqphs3+sqphs3)
c
c.....short channel effect factor
c
200   continue
      if ( (xj.eq.0.0d0).or.(xd.eq.0.0d0) ) go to 210
      wps=xd*sqphbs
      onxj=1.0d0/xj
      xjonxl=xj*onxl
      djonxj=xld*onxj
      wponxj=wps*onxj
      wconxj=coeff0+coeff1*wponxj+coeff2*wponxj*wponxj
      wcs=wconxj*xj
      arga=wconxj+djonxj
      argc=wponxj/(1.0d0+wponxj)
      argb=dsqrt(1.0d0-argc*argc)
      fshort=1.0d0-xjonxl*(arga*argb-djonxj)
      dwpdvb=xd*dsqdvb
      dadvb=(coeff1+coeff2*(wponxj+wponxj))*dwpdvb*onxj
      dbdvb=-argc*argc*(1.0d0-argc)*dwpdvb/(argb*wps)
      dfsdvb=-xjonxl*(dadvb*argb+arga*dbdvb)
      go to 220
210   continue
      fshort=1.0d0
      dfsdvb=0.0d0
      wcs=0.05d-6
c
c.....body effect
c
220   continue
      gammas=gamma*fshort
      fbodys=0.5d0*gammas/(sqphbs+sqphbs)
      fbody=fbodys+fnarrw
      onfbdy=1.0d0/(1.0d0+fbody)
      dfbdvb=-fbodys*dsqdvb/sqphbs+fbodys*dfsdvb/fshort
      qbonco=gammas*sqphbs+fnarrw*phibs
      dqbdvb=gammas*dsqdvb+gamma*dfsdvb*sqphbs-fnarrw
c
c.....static feedback effect
c
      vbix=vbi-eta*vds
c
c.....threshold voltage
c
      vth=vbix+qbonco
      dvtdvd=-eta
      dvtdvb=dqbdvb
c
c.....joint weak inversion and strong inversion
c
      von=vth
      if ( xnfs.eq.0.0d0 ) go to 250
           csonco=charge*xnfs*xl*xw/cox
           cdonco=qbonco/(phibs+phibs)
           xn=1.0d0+csonco+cdonco
           von=vth+vt*xn
           dxndvb=dqbdvb/(phibs+phibs)-qbonco*dsqdvb/(phibs*sqphbs)
           dvodvd=dvtdvd
           dvodvb=dvtdvb+vt*dxndvb
           go to 300
c
c.....cutoff region
c
250   continue
      if ( vgs.gt.von ) go to 300
      cdrain=0.0d0
      gm=0.0d0
      gds=0.0d0
      gmbs=0.0d0
      if ( icharg.ne.0 ) go to 800
      go to 1000
c
c.....device is on
c
300   continue
      vgsx=dmax1(vgs,von)
c
c.....mobility modulation by gate voltage
c
      onfg=1.0d0+theta*(vgsx-vth)
      fgate=1.0d0/onfg
      us=uo*fgate
      dfgdvg=-theta*fgate*fgate
      dfgdvd=-dfgdvg*dvtdvd
      dfgdvb=-dfgdvg*dvtdvb
c
c.....saturation voltage
c
      vdsat=(vgsx-vth)*onfbdy
      vpof=vdsat
      if ( vmax.gt.0.0d0 ) go to 310
      dvsdvg=onfbdy
      dvsdvd=-dvsdvg*dvtdvd
      dvsdvb=-dvsdvg*dvtdvb-vdsat*dfbdvb*onfbdy
      go to 400
  310 vdsc=xl*vmax/us
      onvdsc=1.0d0/vdsc
      arga=(vgsx-vth)*onfbdy
      argb=dsqrt(arga*arga+vdsc*vdsc)
      vdsat=arga+vdsc-argb
      dvsdga=(1.0d0-arga/argb)*onfbdy
      dvsdvg=dvsdga-(1.0d0-vdsc/argb)*vdsc*dfgdvg*onfg
      dvsdvd=-dvsdvg*dvtdvd
      dvsdvb=-dvsdvg*dvtdvb-arga*dvsdga*dfbdvb
c
c.....current factors in linear region
c
400   continue
      vdsx=dmin1(vds,vdsat)
      if ( vdsx.eq.0.0d0 ) go to 900
      cdo=vgsx-vth-0.5d0*(1.0d0+fbody)*vdsx
      dcodvg=1.0d0
      if (vds.lt.vdsat) dcodvd=-dvtdvd-0.5d0*(1.0d0+fbody)
      dcodvb=-dvtdvb-0.5d0*dfbdvb*vdsx
c
c.....normalized drain current
c
410   continue
      cdnorm=cdo*vdsx
      gm=vdsx
      gds=vgsx-vth-(1.0d0+fbody+dvtdvd)*vdsx
      gmbs=dcodvb*vdsx
c
c.....drain current without velocity saturation effect
c
      cd1=beta*cdnorm
      beta=beta*fgate
      cdrain=beta*cdnorm
      gm=beta*gm+dfgdvg*cd1
      gds=beta*gds+dfgdvd*cd1
      gmbs=beta*gmbs
c
c.....velocity saturation factor
c
      if ( vmax.eq.0.0d0 ) go to 500
      fdrain=1.0d0/(1.0d0+vdsx*onvdsc)
      fd2=fdrain*fdrain
      arga=fd2*vdsx*onvdsc*onfg
      dfddvg=-dfgdvg*arga
      dfddvd=-dfgdvd*arga-fd2*onvdsc
      dfddvb=-dfgdvb*arga
c
c.....drain current
c
      gm=fdrain*gm+dfddvg*cdrain
      gds=fdrain*gds+dfddvd*cdrain
      gmbs=fdrain*gmbs+dfddvb*cdrain
      cdrain=fdrain*cdrain
      beta=beta*fdrain
c
c.....channel length modulation
c
500   continue
      if ( vds.le.vdsat ) go to 700
      if ( vmax.eq.0.0d0 ) go to 510
      if (alpha.eq.0.0d0) go to 700
      cdsat=cdrain
      gdsat=cdsat*(1.0d0-fdrain)*onvdsc
      gdsat=dmax1(1.0d-12,gdsat)
      gdoncd=gdsat/cdsat
      gdonfd=gdsat/(1.0d0-fdrain)
      gdonfg=gdsat*onfg
      dgdvg=gdoncd*gm-gdonfd*dfddvg+gdonfg*dfgdvg
      dgdvd=gdoncd*gds-gdonfd*dfddvd+gdonfg*dfgdvd
      dgdvb=gdoncd*gmbs-gdonfd*dfddvb+gdonfg*dfgdvb
c
      emax=cdsat*onxl/gdsat
      emoncd=emax/cdsat
      emongd=emax/gdsat
      demdvg=emoncd*gm-emongd*dgdvg
      demdvd=emoncd*gds-emongd*dgdvd
      demdvb=emoncd*gmbs-emongd*dgdvb
c
      arga=0.5d0*emax*alpha
      argc=xkappa*alpha
      argb=dsqrt(arga*arga+argc*(vds-vdsat))
      delxl=argb-arga
      dldvd=argc/(argb+argb)
      dldem=0.5d0*(arga/argb-1.0d0)*alpha
      ddldvg=dldem*demdvg
      ddldvd=dldem*demdvd-dldvd
      ddldvb=dldem*demdvb
      go to 520
510   continue
      delxl=dsqrt(xkappa*(vds-vdsat)*alpha)
      dldvd=0.5d0*delxl/(vds-vdsat)
      ddldvg=0.0d0
      ddldvd=-dldvd
      ddldvb=0.0d0
c
c.....punch through approximation
c
520   continue
      if ( delxl.le.(0.5d0*xl) ) go to 600
      wcs2=wcs*wcs
      delxl=xl-(xl**2/(4.0d0*delxl))
      arga=4.0d0*(xl-delxl)**2/xl**2
      ddldvg=ddldvg*arga
      ddldvd=ddldvd*arga
      ddldvb=ddldvb*arga
       dldvd= dldvd*arga
c
c.....saturation region
c
600   continue
      dlonxl=delxl*onxl
      xlfact=1.0d0/(1.0d0-dlonxl)
      cdrain=cdrain*xlfact
      diddl=cdrain/(xl-delxl)
      gm=gm*xlfact+diddl*ddldvg
      gds0=gds*xlfact+diddl*ddldvd
      gmbs=gmbs*xlfact+diddl*ddldvb
      gm=gm+gds0*dvsdvg
      gmbs=gmbs+gds0*dvsdvb
      gds=gds0*dvsdvd+diddl*dldvd
c
c.....finish strong inversion case
c
700   continue
      if ( vgs.ge.von ) go to 750
c
c.....weak inversion
c
                onxn=1.0d0/xn
                ondvt=onxn/vt
                wfact=dexp( (vgs-von)*ondvt )
                cdrain=cdrain*wfact
                gms=gm*wfact
                gmw=cdrain*ondvt
                gm=gmw
                if (vds.gt.vdsat) gm=gm+gds0*dvsdvg*wfact
                gds=gds*wfact+(gms-gmw)*dvodvd
                gmbs=gmbs*wfact+(gms-gmw)*dvodvb
     1                         -gmw*(vgs-von)*onxn*dxndvb
c
c.....charge computation
c
  750 continue
      if (icharg.eq.0) go to 1000
      if (vgs.le.vth) go to 800
      call mqspof(vds,vbs,vgs,vpof,vdsat1,vth,vbin,gamasd,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      go to 2000
c
c.....charge computation for vgs<vth
c
800   continue
      xqc=xqco
      call mosq3(vds,vbs,vpof,vdsat1,vth,vbin,gamasd,cox,phi,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      qspof=0.0d0
      go to 2000
c
c.....special case of vds=0.0d0
c
900   continue
      beta=beta*fgate
      cdrain=0.0d0
      gm=0.0d0
      gds=beta*(vgsx-vth)
      gmbs=0.0d0
           if ( (xnfs.ne.0.0d0).and.(vgs.lt.von) )
     1          gds=gds*dexp((vgs-von)/(vt*xn))
      if (icharg.eq.0) go to 1000
      call mosq3(vds,vbs,vpof,vdsat1,vth,vbin,gamasd,cox,phi,
     1   qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
1000  qspof=0.0d0
c
c.....done
c
 2000 return
      end
