      subroutine nlcsrc
      implicit double precision (a-h,o-z)
c
c     this routine loads the nonlinear controlled sources into the
c coefficient matrix.
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
c  nonlinear voltage-controlled current sources
c
      loc=locate(5)
   10 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 100
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      lnod=nodplc(loc+6)
      lmat=nodplc(loc+7)
      lcoef=nodplc(loc+8)
      call sizmem(nodplc(loc+8),ncoef)
      larg=nodplc(loc+9)
      lexp=nodplc(loc+10)
      lic=nodplc(loc+11)
      loct=nodplc(loc+12)+1
      icheck=0
      do 20 i=1,ndim
      call update(value(lic+i),loct,nodplc(lnod+1),nodplc(lnod+2),2,
     1   icheck)
      value(larg+i)=value(lx0+loct)
      loct=loct+2
      lnod=lnod+2
   20 continue
      call evpoly(cold,0,lcoef,ncoef,larg,ndim,lexp)
      loct=nodplc(loc+12)
      if (icheck.eq.1) go to 30
      if (initf.eq.6) go to 30
      tol=reltol*dmax1(dabs(cold),dabs(value(lx0+loct)))+abstol
      if (dabs(cold-value(lx0+loct)).lt.tol) go to 40
   30 noncon=noncon+1
   40 value(lx0+loct)=cold
      ceq=cold
      do 50 i=1,ndim
      call evpoly(geq,i,lcoef,ncoef,larg,ndim,lexp)
      loct=loct+2
      value(lx0+loct)=geq
      ceq=ceq-geq*value(larg+i)
      locy=lvn+nodplc(lmat+1)
      value(locy)=value(locy)+geq
      locy=lvn+nodplc(lmat+2)
      value(locy)=value(locy)-geq
      locy=lvn+nodplc(lmat+3)
      value(locy)=value(locy)-geq
      locy=lvn+nodplc(lmat+4)
      value(locy)=value(locy)+geq
      lmat=lmat+4
   50 continue
      value(lvn+node1)=value(lvn+node1)-ceq
      value(lvn+node2)=value(lvn+node2)+ceq
      loc=nodplc(loc)
      go to 10
c
c  nonlinear voltage controlled voltage sources
c
  100 loc=locate(6)
  110 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 200
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      iptr=nodplc(loc+6)
      lnod=nodplc(loc+7)
      lmat=nodplc(loc+8)
      lcoef=nodplc(loc+9)
      call sizmem(nodplc(loc+9),ncoef)
      larg=nodplc(loc+10)
      lexp=nodplc(loc+11)
      lic=nodplc(loc+12)
      loct=nodplc(loc+13)+2
      icheck=0
      do 120 i=1,ndim
      call update(value(lic+i),loct,nodplc(lnod+1),nodplc(lnod+2),2,
     1   icheck)
      value(larg+i)=value(lx0+loct)
      loct=loct+2
      lnod=lnod+2
  120 continue
      call evpoly(volt,0,lcoef,ncoef,larg,ndim,lexp)
      loct=nodplc(loc+13)
      if (icheck.eq.1) go to 130
      if (initf.eq.6) go to 130
      tol=reltol*dmax1(dabs(volt),dabs(value(lx0+loct)))+vntol
      if (dabs(volt-value(lx0+loct)).lt.tol) go to 140
  130 noncon=noncon+1
  140 value(lx0+loct)=volt
      value(lx0+loct+1)=value(lvnim1+iptr)
      veq=volt
      locy=lvn+nodplc(lmat+1)
      value(locy)=+1.0d0
      locy=lvn+nodplc(lmat+2)
      value(locy)=-1.0d0
      locy=lvn+nodplc(lmat+3)
      value(locy)=+1.0d0
      locy=lvn+nodplc(lmat+4)
      value(locy)=-1.0d0
      lmat=lmat+4
      loct=loct+1
      do 150 i=1,ndim
      call evpoly(vgain,i,lcoef,ncoef,larg,ndim,lexp)
      loct=loct+2
      value(lx0+loct)=vgain
      veq=veq-vgain*value(larg+i)
      locy=lvn+nodplc(lmat+1)
      value(locy)=value(locy)-vgain
      locy=lvn+nodplc(lmat+2)
      value(locy)=value(locy)+vgain
      lmat=lmat+2
  150 continue
      value(lvn+iptr)=veq
      loc=nodplc(loc)
      go to 110
c
c  nonlinear current-controlled current sources
c
  200 loc=locate(7)
  210 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 300
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      lvs=nodplc(loc+6)
      lmat=nodplc(loc+7)
      lcoef=nodplc(loc+8)
      call sizmem(nodplc(loc+8),ncoef)
      larg=nodplc(loc+9)
      lexp=nodplc(loc+10)
      lic=nodplc(loc+11)
      loct=nodplc(loc+12)+1
      icheck=0
      do 220 i=1,ndim
      iptr=nodplc(lvs+i)
      iptr=nodplc(iptr+6)
      call update(value(lic+i),loct,iptr,1,2,icheck)
      value(larg+i)=value(lx0+loct)
      loct=loct+2
  220 continue
      call evpoly(csrc,0,lcoef,ncoef,larg,ndim,lexp)
      loct=nodplc(loc+12)
      if (icheck.eq.1) go to 230
      if (initf.eq.6) go to 230
      tol=reltol*dmax1(dabs(csrc),dabs(value(lx0+loct)))+abstol
      if (dabs(csrc-value(lx0+loct)).lt.tol) go to 240
  230 noncon=noncon+1
  240 value(lx0+loct)=csrc
      ceq=csrc
      do 250 i=1,ndim
      call evpoly(cgain,i,lcoef,ncoef,larg,ndim,lexp)
      loct=loct+2
      value(lx0+loct)=cgain
      ceq=ceq-cgain*value(larg+i)
      locy=lvn+nodplc(lmat+1)
      value(locy)=value(locy)+cgain
      locy=lvn+nodplc(lmat+2)
      value(locy)=value(locy)-cgain
      lmat=lmat+2
  250 continue
      value(lvn+node1)=value(lvn+node1)-ceq
      value(lvn+node2)=value(lvn+node2)+ceq
      loc=nodplc(loc)
      go to 210
c
c  nonlinear current controlled voltage sources
c
  300 loc=locate(8)
  310 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 1000
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      ibr=nodplc(loc+6)
      lvs=nodplc(loc+7)
      lmat=nodplc(loc+8)
      lcoef=nodplc(loc+9)
      call sizmem(nodplc(loc+9),ncoef)
      larg=nodplc(loc+10)
      lexp=nodplc(loc+11)
      lic=nodplc(loc+12)
      loct=nodplc(loc+13)+2
      icheck=0
      do 320 i=1,ndim
      iptr=nodplc(lvs+i)
      iptr=nodplc(iptr+6)
      call update(value(lic+i),loct,iptr,1,2,icheck)
      value(larg+i)=value(lx0+loct)
      loct=loct+2
  320 continue
      call evpoly(volt,0,lcoef,ncoef,larg,ndim,lexp)
      loct=nodplc(loc+13)
      if (icheck.eq.1) go to 330
      if (initf.eq.6) go to 330
      tol=reltol*dmax1(dabs(volt),dabs(value(lx0+loct)))+vntol
      if (dabs(volt-value(lx0+loct)).lt.tol) go to 340
  330 noncon=noncon+1
  340 value(lx0+loct)=volt
      value(lx0+loct+1)=value(lvnim1+ibr)
      veq=volt
      locy=lvn+nodplc(lmat+1)
      value(locy)=+1.0d0
      locy=lvn+nodplc(lmat+2)
      value(locy)=-1.0d0
      locy=lvn+nodplc(lmat+3)
      value(locy)=+1.0d0
      locy=lvn+nodplc(lmat+4)
      value(locy)=-1.0d0
      lmat=lmat+4
      loct=loct+1
      do 350 i=1,ndim
      call evpoly(transr,i,lcoef,ncoef,larg,ndim,lexp)
      loct=loct+2
      value(lx0+loct)=transr
      veq=veq-transr*value(larg+i)
      locy=lvn+nodplc(lmat+i)
      value(locy)=value(locy)-transr
  350 continue
      value(lvn+ibr)=veq
      loc=nodplc(loc)
      go to 310
c
c  finished
c
 1000 return
      end
