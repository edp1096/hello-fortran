      subroutine intgr8(geq,ceq,capval,loct)
      implicit double precision (a-h,o-z)
c
c     this routine performs the actual numerical integration for each
c circuit element.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension qcap(1),ccap(1)
      equivalence (qcap(1),value(1)),(ccap(1),value(2))
c
c
      if (method.eq.2) go to 100
c
c  trapezoidal algorithm
c
      if (iord.eq.1) go to 100
      ccap(lx0+loct)=-ccap(lx1+loct)*ag(2)
     1   +ag(1)*(qcap(lx0+loct)-qcap(lx1+loct))
      go to 190
c
c  gears algorithm
c
  100 go to (110,120,130,140,150,160), iord
  110 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
      go to 190
  120 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
     1              +ag(3)*qcap(lx2+loct)
      go to 190
  130 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
     1              +ag(3)*qcap(lx2+loct)+ag(4)*qcap(lx3+loct)
      go to 190
  140 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
     1              +ag(3)*qcap(lx2+loct)+ag(4)*qcap(lx3+loct)
     2              +ag(5)*qcap(lx4+loct)
      go to 190
  150 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
     1              +ag(3)*qcap(lx2+loct)+ag(4)*qcap(lx3+loct)
     2              +ag(5)*qcap(lx4+loct)+ag(6)*qcap(lx5+loct)
      go to 190
  160 ccap(lx0+loct)=ag(1)*qcap(lx0+loct)+ag(2)*qcap(lx1+loct)
     1              +ag(3)*qcap(lx2+loct)+ag(4)*qcap(lx3+loct)
     2              +ag(5)*qcap(lx4+loct)+ag(6)*qcap(lx5+loct)
     3              +ag(7)*qcap(lx6+loct)
c... ceq is the equivalent current applicable to linear capacitance
c    (inductance) only, i.e. q=c*v
  190 ceq=ccap(lx0+loct)-ag(1)*qcap(lx0+loct)
      geq=ag(1)*capval
      return
      end
