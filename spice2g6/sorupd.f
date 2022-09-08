      subroutine sorupd
      implicit double precision (a-h,o-z)
c
c     this routine updates the independent voltage and current sources
c used in the circuit.  it also updates the ltd table (which contains
c previous (delayed) values of the sources used to model transmission
c lines).
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      do 500 id=9,10
      loc=locate(id)
   10 if (loc.eq.0) go to 500
      if ((id.eq.9).and.(nodplc(loc+11).ne.0)) go to 500
      if ((id.eq.10).and.(nodplc(loc+6).ne.0)) go to 500
      locv=nodplc(loc+1)
      locp=nodplc(loc+5)
      itype=nodplc(loc+4)+1
      go to (490,100,200,300,400,450), itype
c
c  pulse source
c
  100 v1=value(locp+1)
      v2=value(locp+2)
      t1=value(locp+3)
      t2=value(locp+4)
      t3=value(locp+5)
      t4=value(locp+6)
      period=value(locp+7)
      time1=time
      if (time1.le.0.0d0) go to 160
  110 if (time1.lt.t1+period) go to 120
      time1=time1-period
      go to 110
  120 if (time1.lt.t4) go to 130
      value(locv+1)=v1
      go to 490
  130 if (time1.lt.t3) go to 140
      value(locv+1)=v2+(time1-t3)*(v1-v2)/(t4-t3)
      go to 490
  140 if (time1.lt.t2) go to 150
      value(locv+1)=v2
      go to 490
  150 if (time1.lt.t1) go to 160
      value(locv+1)=v1+(time1-t1)*(v2-v1)/(t2-t1)
      go to 490
  160 value(locv+1)=v1
      go to 490
c
c  sinusoidal source
c
  200 v1=value(locp+1)
      v2=value(locp+2)
      omeg=value(locp+3)
      t1=value(locp+4)
      theta=value(locp+5)
      time1=time-t1
      if (time1.gt.0.0d0) go to 210
      value(locv+1)=v1
      go to 490
  210 if (theta.ne.0.0d0) go to 220
      value(locv+1)=v1+v2*dsin(omeg*time1)
      go to 490
  220 value(locv+1)=v1+v2*dsin(omeg*time1)*dexp(-time1*theta)
      go to 490
c
c  exponential source
c
  300 v1=value(locp+1)
      v2=value(locp+2)
      t1=value(locp+3)
      tau1=value(locp+4)
      t2=value(locp+5)
      tau2=value(locp+6)
      time1=time
      if (time1.gt.t1) go to 310
      value(locv+1)=v1
      go to 490
  310 if (time1.gt.t2) go to 320
      value(locv+1)=v1+(v2-v1)*(1.0d0-dexp((t1-time1)/tau1))
      go to 490
  320 value(locv+1)=v1+(v2-v1)*(1.0d0-dexp((t1-time1)/tau1))
     1   +(v1-v2)*(1.0d0-dexp((t2-time1)/tau2))
      go to 490
c
c  piecewise-linear source
c
  400 t1=value(locp+1)
      v1=value(locp+2)
      t2=value(locp+3)
      v2=value(locp+4)
      iknt=4
  410 if (time.le.t2) go to 420
      t1=t2
      v1=v2
      t2=value(locp+iknt+1)
      v2=value(locp+iknt+2)
      iknt=iknt+2
      go to 410
  420 value(locv+1)=v1+((time-t1)/(t2-t1))*(v2-v1)
      go to 490
c
c  single-frequency fm
c
  450 v1=value(locp+1)
      v2=value(locp+2)
      omegc=value(locp+3)
      xmod=value(locp+4)
      omegs=value(locp+5)
      value(locv+1)=v1+v2*dsin(omegc*time+xmod*dsin(omegs*time))
  490 loc=nodplc(loc)
      go to 10
  500 continue
c
c  update transmission line sources
c
      if (jelcnt(17).eq.0) go to 1000
      if (mode.ne.2) go to 1000
      call sizmem(ltd,ltdsiz)
      numtd=ltdsiz/ntlin
      if (numtd.lt.3) go to 900
      loc=locate(17)
  610 if (loc.eq.0) go to 1000
      locv=nodplc(loc+1)
      td=value(locv+2)
      baktim=time-td
      if (baktim.lt.0.0d0) go to 640
      ltdptr=nodplc(loc+30)
      icntr=2
      l1=ltd
      l2=l1+ntlin
      l3=l2+ntlin
      t1=value(l1+1)
      t2=value(l2+1)
  620 t3=value(l3+1)
      icntr=icntr+1
      if (baktim.le.t3) go to 630
      if (icntr.eq.numtd) go to 900
      l1=l2
      l2=l3
      l3=l2+ntlin
      t1=t2
      t2=t3
      go to 620
  630 dt1t2=t1-t2
      dt1t3=t1-t3
      dt2t3=t2-t3
      tdnom1=1.0d0/(dt1t2*dt1t3)
      tdnom2=-1.0d0/(dt1t2*dt2t3)
      tdnom3=1.0d0/(dt2t3*dt1t3)
      dtt1=baktim-t1
      dtt2=baktim-t2
      dtt3=baktim-t3
      tfact1=dtt2*dtt3*tdnom1
      tfact2=dtt1*dtt3*tdnom2
      tfact3=dtt1*dtt2*tdnom3
      value(locv+3)=value(l1+ltdptr+0)*tfact1+value(l2+ltdptr+0)*tfact2
     1   +value(l3+ltdptr+0)*tfact3
      value(locv+4)=value(l1+ltdptr+1)*tfact1+value(l2+ltdptr+1)*tfact2
     1   +value(l3+ltdptr+1)*tfact3
  640 loc=nodplc(loc)
      go to 610
c
c  internal logic error:  less than 3 entries in ltd
c
  900 nogo=1
      write (iofile,901) numtd,icntr
  901 format('0*abort*:  internal spice error:  sorupd:  ',2i5/)
c
c  finished
c
 1000 return
      end
