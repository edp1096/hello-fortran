      subroutine moscap(vgd,vgs,vgb,covlgd,covlgs,covlgb,
     1   capbd,capbs,cggb,cgdb,cgsb,cbgb,cbdb,cbsb,
     2   gcggb,gcgdb,gcgsb,gcbgb,gcbdb,gcbsb,
     3   gcdgb,gcddb,gcdsb,gcsgb,gcsdb,gcssb,
     4   qgate,qchan,qbulk,qdrn,qsrc)
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c     compute equivalent conductances
c     divide up the channel charge (1-xqc)/xqc to source and drain
c
      gcg=(cggb+cbgb)*ag(1)
      gcd=(cgdb+cbdb)*ag(1)
      gcs=(cgsb+cbsb)*ag(1)
      gcgxd=-xqc*gcg
      gcgxs=-(1.0d0-xqc)*gcg
      gcdxd=-xqc*gcd
      gcdxs=-(1.0d0-xqc)*gcd
      gcsxd=-xqc*gcs
      gcsxs=-(1.0d0-xqc)*gcs
      gcdgb=gcgxd-covlgd*ag(1)
      gcddb=gcdxd+(capbd+covlgd)*ag(1)
      gcdsb=gcsxd
      gcsgb=gcgxs-covlgs*ag(1)
      gcsdb=gcdxs
      gcssb=gcsxs+(capbs+covlgs)*ag(1)
      gcggb=(cggb+covlgd+covlgs+covlgb)*ag(1)
      gcgdb=(cgdb-covlgd)*ag(1)
      gcgsb=(cgsb-covlgs)*ag(1)
      gcbgb=(cbgb-covlgb)*ag(1)
      gcbdb=(cbdb-capbd)*ag(1)
      gcbsb=(cbsb-capbs)*ag(1)
c
c     compute total terminal charges
c
      qgd=covlgd*vgd
      qgs=covlgs*vgs
      qgb=covlgb*vgb
      qgate=qgate+qgd+qgs+qgb
      qbulk=qbulk-qgb
      qdrn=xqc*qchan-qgd
      qsrc=(1.0d0-xqc)*qchan-qgs
c
c     finished
c
      return
      end
