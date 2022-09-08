      subroutine terr(loct,delnew)
      implicit double precision (a-h,o-z)
c
c     this routine estimates the local truncation error for a particular
c circuit element.  it then computes the appropriate stepsize which
c should be used.
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
      dimension qcap(1),ccap(1),diff(8),deltmp(7),coef(6)
      equivalence (qcap(1),value(1)),(ccap(1),value(2))
      data coef / 5.000000000d-1, 2.222222222d-1, 1.363636364d-1,
     1            9.600000000d-2, 7.299270073d-2, 5.830903790d-2 /
      data xtwelv / 8.333333333d-2 /
c
c
      tol=reltol*dmax1(dabs(ccap(lx0+loct)),dabs(ccap(lx1+loct)))+abstol
      ctol=reltol*dmax1(dabs(qcap(lx0+loct)),dabs(qcap(lx1+loct)),
     1   chgtol)/delta
      tol=dmax1(tol,ctol)
c
c  determine divided differences
c
      go to (6,5,4,3,2,1), iord
    1 diff(8)=qcap(lx7+loct)
    2 diff(7)=qcap(lx6+loct)
    3 diff(6)=qcap(lx5+loct)
    4 diff(5)=qcap(lx4+loct)
    5 diff(4)=qcap(lx3+loct)
    6 diff(3)=qcap(lx2+loct)
      diff(2)=qcap(lx1+loct)
      diff(1)=qcap(lx0+loct)
      istop=iord+1
      do 10 i=1,istop
      deltmp(i)=delold(i)
   10 continue
   20 do 30 i=1,istop
      diff(i)=(diff(i)-diff(i+1))/deltmp(i)
   30 continue
      istop=istop-1
      if (istop.eq.0) go to 100
      do 40 i=1,istop
      deltmp(i)=deltmp(i+1)+delold(i)
   40 continue
      go to 20
c
c  diff(1) contains divided difference
c
  100 const=coef(iord)
      if ((method.eq.1).and.(iord.eq.2)) const=xtwelv
      del=trtol*tol/dmax1(abstol,const*dabs(diff(1)))
      if (iord.eq.1) go to 200
      if (iord.ge.3) go to 150
      del=dsqrt(del)
      go to 200
  150 del=dexp(dlog(del)/dble(iord))
  200 delnew=dmin1(delnew,del)
      return
      end
