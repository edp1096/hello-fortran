      subroutine reserv (node1,node2)
      implicit double precision (a-h,o-z)
c
c     this routine records the fact that the (node1, node2) element of
c the circuit equation coefficient matrix is nonzero.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      logical memptr
c
      if (nogo.ne.0) go to 300
c...  test for ground
      if (node1.eq.1) go to 300
      if (node2.eq.1) go to 300
c
c     reserve (node1,node2) in row node1 at col posn node2
c
      loc=node1
   10 locj=loc
      loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 20
      if (nodplc(jcolno+loc)-node2) 10,300,20
   20 call sizmem(jcpt,isize)
      newloc=isize+1
      nodplc(numoff+node1)=nodplc(numoff+node1)+1
      nodplc(nmoffc+node2)=nodplc(nmoffc+node2)+1
      call extmem(jcpt,1)
      call extmem(jcolno,1)
      nodplc(jcpt+locj)=newloc
      nodplc(jcpt+newloc)=loc
      nodplc(jcolno+newloc)=node2
c
c     reserve (node1,node2) in col node2 at row posn node1
c
      loc=node2
   30 loci=loc
      loc=nodplc(irpt+loc)
      if (loc.eq.0) go to 40
      if (nodplc(irowno+loc)-node1) 30,300,40
   40 call extmem(irpt,1)
      call extmem(irowno,1)
      nodplc(irpt+loci)=newloc
      nodplc(irpt+newloc)=loc
      nodplc(irowno+newloc)=node1
c
c     mark diagonal
c
      if (node1.ne.node2) go to 300
      if (memptr(ndiag)) nodplc(ndiag+node1)=1
c
c     finished
c
  300 return
      end
