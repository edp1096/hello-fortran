      subroutine asol
      implicit double precision (a-h,o-z)
c
c     this routine evaluates the adjoint circuit response by doing a
c forward/backward substitution on the transpose of the coefficient
c matrix.
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
c  forward substitution
c
      do 20 i=2,nstop
      iord=nodplc(icswpf+i)
      loc=i
   10 loc=nodplc(irpt+loc)
      if (nodplc(irowno+loc).ge.i) go to 15
      j=nodplc(irowno+loc)
      jord=nodplc(icswpf+j)
      value(lvn+iord)=value(lvn+iord)-value(lvn+loc)*value(lvn+jord)
      go to 10
   15 jord=nodplc(irswpf+i)
      locnn=indxx(jord,iord)
      value(lvn+iord)=value(lvn+iord)/value(lvn+locnn)
   20 continue
c
c  backward substitution
c
      i=nstop
   30 i=i-1
      if (i.le.1) go to 60
      iord=nodplc(icswpf+i)
      loc=i
   35 loc=nodplc(irpt+loc)
   40 if (nodplc(irowno+loc).ne.i) go to 35
   50 loc=nodplc(irpt+loc)
      if (loc.eq.0) go to 30
      j=nodplc(irowno+loc)
      jord=nodplc(icswpf+j)
      value(lvn+iord)=value(lvn+iord)-value(lvn+loc)*value(lvn+jord)
      go to 50
c
c     reorder solution vector
c
   60 do 70 i=1,nstop
      j=nodplc(irswpr+i)
      k=nodplc(icswpf+j)
      value(lvntmp+i)=value(lvn+k)
   70 continue
      call copy8(value(lvntmp+1),value(lvn+1),nstop)
c
c  finished
c
      return
      end
