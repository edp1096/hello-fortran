      subroutine mqspof(vds,vbs,vgs,vpof,vdsat1,vth,vbin,gamasd,
     $qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      implicit double precision (a-h,o-z)
c
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c
c     vdsat1=dmax1(vds,vdsat1)+1.0d-3
      if( lev .eq. 3 ) goto 50
      if( lev .ne. 2 ) goto 1000
      call mosq2(vds,vbs,vgs,vdsat,vth,vbin,gamasd,cox,phi,
     $qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
      if (vds.ge.vdsat) go to 80
      call mosq2(vds,vbs,vpof,vdsat1,vth,vbin, gamasd,cox,phi,
     $qg1,qcpof1,qb1,cggb1,cgdb1,cgsb1,cbgb1,cbdb1,cbsb1)
      call mosq2(vdsat,vbs,vgs,vdsat,vth,vbin, gamasd,cox,phi,
     $qg2,qcpof2,qb2,cggb2,cgdb2,cgsb2,cbgb2,cbdb2,cbsb2)
      goto 75
   50 call mosq3(vds,vbs,vpof,vdsat1,vth,vbin, gamasd,cox,phi,
     $qg,qcpof,qb,cggb1,cgdb1,cgsb1,cbgb1,cbdb1,cbsb1)
      call mosq3(vds,vbs,vgs,vdsat,vth,vbin,gamasd,cox,phi,
     $qg,qc,qb,cggb,cgdb,cgsb,cbgb,cbdb,cbsb)
   75 if(vgs.gt.vpof. or .vds.lt.vdsat) goto 100
   80 xqc = xqco
      goto 1000
c
c     tangential limiting of qs
c
  100 csgb1=-(1.0d0-xqco)*(cggb1+cbgb1)
      qs=csgb1*(vgs-vpof)
     1   +(1.0d0-xqco)*qcpof1
c      write(iofile,*) "vgs,vds,qc,cggb,cgdb,cgsb,cbgb,cbdb,cbsb =",
c    1   vgs,vds,qc,cggb,cgdb,cgsb,cbgb,cbdb,cbsb
c      write(iofile,*) "vpof,vdsat,vdsat1,qcpof1,qcpof2,qs,csgb1 =",
c    1   vpof,vdsat,vdsat1,qcpof1,qcpof2,qs,csgb1
      qspof2=(1.0d0-xqco)*qcpof2
      if (dabs(qs) .lt. dabs(qspof2)) qs=qspof2
      if( dabs( qs ) .ge. 0.5d0 * dabs( qc ) ) goto 200
c     csdb=-0.25d0*(cgdb+cbdb)
c     qs=qs+csdb*(vdsat-vds)
c     xqc=dmin1(0.5d0,(qc-qs)/qc)
      xqc=0.5d0
c      write(iofile,*) "qs,xqc =",
c    1   qs,xqc
      goto 1000
  200 qd = qc - qs
      xqc = qd / qc
c     write(iofile,*) "200,qs,qd,xqc =",
c    1   qs,qd,xqc
c
c     constant limiting of qs
c
c 100 qdpof = qcpof * xqco
c     qspof = qcpof - qdpof
c     if( dabs( qspof ) .gt. 0.5d0 * dabs( qc ) ) goto 200
c     xqc = 0.5d0
c     goto 1000
c 200 qd = qc - qspof
c     qs = qspof
c     xqc = qd / qc
 1000 return
      end
