      integer function indxx(node1,node2)
      implicit double precision (a-h,o-z)
c
c     this routine maps a (row, column) matrix term specification into
c the offset from the origin of the matrix storage at which the term is
c actually located.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      if (node1.eq.1) go to 100
      if (node2.eq.1) go to 100
c
      n1=nodplc(irswpr+node1)
      n2=nodplc(icswpr+node2)
      if (n1-n2) 10,10,30
c
c     search col n2
c
   10 loc=n2
   15 loc=nodplc(irpt+loc)
      if (loc.eq.0) go to 100
      if (nodplc(irowno+loc)-n1) 15,20,15
   20 indxx=loc
      return
c
c     search row n1
c
   30 loc=n1
   35 loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 100
      if (nodplc(jcolno+loc)-n2) 35,40,35
   40 indxx=loc
      return
c
c     unused location
c
  100 indxx=1
      return
      end
