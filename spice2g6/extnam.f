      subroutine extnam(aname,index)
      implicit double precision (a-h,o-z)
c
c     this routine adds 'aname' to the list of 'unsatisfied' names (that
c is, names which can only be resolved after subcircuit expansion).
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      integer xxor
c
c
      anam=aname
      if (nunsat.eq.0) go to 20
      do 10 index=1,nunsat
      if (xxor(anam,value(iunsat+index)).eq.0) go to 30
   10 continue
c
   20 call extmem(iunsat,1)
      nunsat=nunsat+1
      index=nunsat
      value(iunsat+index)=anam
   30 return
      end
