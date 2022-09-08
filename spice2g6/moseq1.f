      subroutine moseq1(vds,vbs,vgs,gm,gds,gmbs)
      implicit double precision (a-h,o-z)
c
c     this routine evaluates the drain current and its derivatives
c     using the shichman-hodges model and the charges associated
c     with the gate, channel and bulk for mosfets
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
      vbd=vbs-vds
      vgb=vgs-vbs
c
c
      if (vbs.gt.0.0d0) go to 102
      sarg=dsqrt(phi-vbs)
      go to 104
  102 sarg=dsqrt(phi)
      sarg=sarg-vbs/(sarg+sarg)
      sarg=dmax1(0.0d0,sarg)
  104 von=vbi+gamma*sarg
      vgst=vgs-von
      vdsat=dmax1(vgst,0.0d0)
      if (sarg.gt.0.0d0) go to 105
      arg=0.0d0
      go to 108
  105 arg=gamma/(sarg+sarg)
  108 if (vgst.gt.0.0d0) go to 110
c
c     cutoff region
c
      cdrain=0.0d0
      gm=0.0d0
      gds=0.0d0
      gmbs=0.0d0
      go to 1000
c
c     saturation region
c
  110 betap=beta*(1.0d0+xlamda*vds)
      if (vgst.gt.vds) go to 120
      cdrain=betap*vgst*vgst*0.5d0
      gm=betap*vgst
      gds=xlamda*beta*vgst*vgst*0.5d0
      gmbs=gm*arg
      go to 1000
c
c     linear region
c
  120 cdrain=betap*vds*(vgst-0.5d0*vds)
      gm=betap*vds
      gds=betap*(vgst-vds)+xlamda*beta*vds*(vgst-0.5d0*vds)
      gmbs=gm*arg
c
c     finished
c
 1000 return
      end
