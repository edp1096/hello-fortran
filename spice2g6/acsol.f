      subroutine acsol
      implicit double precision (a-h,o-z)
c
c     this routine solves the circuit equations by performing a forward
c and backward substitution using the previously-computed lu factors.
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
      loc=i
      iord=nodplc(irswpf+i)
   10 loc=nodplc(jcpt+loc)
      if (nodplc(jcolno+loc).ge.i) go to 20
      j=nodplc(jcolno+loc)
      jord=nodplc(irswpf+j)
      call cmult(value(lynl+loc),value(imynl+loc),
     1     value(lvn+jord),value(imvn+jord),xreal,ximag)
      value(lvn+iord)=value(lvn+iord)-xreal
      value(imvn+iord)=value(imvn+iord)-ximag
      go to 10
   20 continue
c
c      back substitution
c
      i=nstop
      iord=nodplc(irswpf+i)
      jord=nodplc(icswpf+i)
      locnn=indxx(iord,jord)
   30 call cdiv(value(lvn+iord),value(imvn+iord),value(lynl+locnn),
     1     value(imynl+locnn),value(lvn+iord),value(imvn+iord))
      i=i-1
      if (i.le.1) go to 60
      iord=nodplc(irswpf+i)
      loc=i
   35 loc=nodplc(jcpt+loc)
   40 if (nodplc(jcolno+loc).ne.i) go to 35
      locnn=loc
   50 loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 30
      j=nodplc(jcolno+loc)
      jord=nodplc(irswpf+j)
      call cmult(value(lynl+loc),value(imynl+loc),
     1     value(lvn+jord),value(imvn+jord),xreal,ximag)
      value(lvn+iord)=value(lvn+iord)-xreal
      value(imvn+iord)=value(imvn+iord)-ximag
      go to 50
c
c  reorder solution vector
c
   60 do 70 i=1,nstop
      j=nodplc(icswpr+i)
      k=nodplc(irswpf+j)
      value(ndiag+i)=value(lvn+k)
      value(ndiag+i+nstop)=value(imvn+k)
   70 continue
      call copy8(value(ndiag+1),value(lvn+1),nstop)
      call copy8(value(ndiag+1+nstop),value(imvn+1),nstop)
      do 120 i=2,nstop
      cvalue(lcvn+i)=cmplx(sngl(value(lvn+i)),sngl(value(imvn+i)))
  120 continue
      cvalue(lcvn+1)=cmplx(0.0e0,0.0e0)
c
c  finished
c
      return
      end
