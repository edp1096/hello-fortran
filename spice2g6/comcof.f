      subroutine comcof
      implicit double precision (a-h,o-z)
c
c     this routine calculates the timestep-dependent terms used in the
c numerical integration.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      dimension gmat(7,7)
c
c  compute coefficients for particular integration method
c
      if (method.ne.1) go to 5
      if (iord.eq.1) go to 5
c...  trapezoidal method
      ag(1)=1.0d0/delta/(1.0d0-xmu)
      ag(2)=xmu/(1.0d0-xmu)
      go to 200
c
c  construct gear coefficient matrix
c
    5 istop=iord+1
      call zero8(ag,istop)
      ag(2)=-1.0d0
      do 10 i=1,istop
      gmat(1,i)=1.0d0
   10 continue
      do 20 i=2,istop
      gmat(i,1)=0.0d0
   20 continue
      arg=0.0d0
      do 40 i=2,istop
      arg=arg+delold(i-1)
      arg1=1.0d0
      do 30 j=2,istop
      arg1=arg1*arg
      gmat(j,i)=arg1
   30 continue
   40 continue
c
c  solve for gear coefficients ag(*)
c
c
c  lu decomposition
c
      do 70 i=2,istop
      jstart=i+1
      if (jstart.gt.istop) go to 70
      do 60 j=jstart,istop
      gmat(j,i)=gmat(j,i)/gmat(i,i)
      do 50 k=jstart,istop
      gmat(j,k)=gmat(j,k)-gmat(j,i)*gmat(i,k)
   50 continue
   60 continue
   70 continue
c
c  forward substitution
c
      do 90 i=2,istop
      jstart=i+1
      if (jstart.gt.istop) go to 90
      do 80 j=jstart,istop
      ag(j)=ag(j)-gmat(j,i)*ag(i)
   80 continue
   90 continue
c
c  backward substitution
c
      ag(istop)=ag(istop)/gmat(istop,istop)
      ir=istop
      do 110 i=2,istop
      jstart=ir
      ir=ir-1
      do 100 j=jstart,istop
      ag(ir)=ag(ir)-gmat(ir,j)*ag(j)
  100 continue
      ag(ir)=ag(ir)/gmat(ir,ir)
  110 continue
c
c  finished
c
  200 return
      end
