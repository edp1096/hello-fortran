      subroutine cmeyer (vgs0,vgd0,vgb0,von0,vdsat0,vgs1,vgd1,vgb1,
     1   covlgs,covlgd,covlgb,cgs0,cgd0,cgb0,cgs1,cgd1,cgb1)
      implicit double precision (a-h,o-z)
c
c     this routine computes the mosfet overlap capacitances as functions
c of the device terminal voltages.
c
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
      indax=1
      vgs=vgs1
      vgd=vgd1
      vgb=vgb1
      vons=von
      vbs=vgs-vgb
      vdbsat=vdsat-vbs
      vdb=vgb-vgd
   10 vds=vgs-vgd
      vgbt=vgs-vons
      if (vgbt.gt.-phi) go to 100
      cgb=cox+covlgb
      cgd=covlgd
      cgs=covlgs
      go to 430
c
c
  100 if (vgbt.gt.-phi/2.0d0) go to 200
      cgb=-vgbt*cox/phi+covlgb
      cgd=covlgd
      cgs=covlgs
      go to 430
c
c
  200 if (vgbt.gt.0.0d0) go to 300
      cgb=-vgbt*cox/phi+covlgb
      cgd=covlgd
      cgs=cox/(7.5d-1*phi)*vgbt+cox/1.5d0+covlgs
      go to 430
c
c
  300 if (vdbsat.gt.vdb) go to 400
      cgb=covlgb
      cgd=covlgd
      cgs=cox/1.5d0+covlgs
      go to 430
c
c
  400 vddif=2.0d0*vdbsat-vdb
      vddif1=vdbsat-vdb-1.0d-12
      vddif2=vddif*vddif
      cgd=cox*(1.0d0-vdbsat*vdbsat/vddif2)/1.5d0+covlgd
      cgs=cox*(1.0d0-vddif1*vddif1/vddif2)/1.5d0+covlgs
      cgb=covlgb
c
c
  430 go to (440,560), indax
  440 indax=2
      cgs1=cgs
      cgd1=cgd
      cgb1=cgb
      vgs=vgs0
      vgd=vgd0
      vgb=vgb0
      vons=von0
      vbs=vgs-vgb
      vdbsat=vdsat0-vbs
      vdb=vgb-vgd
      go to 10
c
c
  560 cgs0=cgs
      cgd0=cgd
      cgb0=cgb
c
c  finished
c
 1000 return
      end
