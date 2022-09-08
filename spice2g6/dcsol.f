      subroutine dcsol
      implicit double precision (a-h,o-z)
c
c     this routine solves the system of circuit equations by performing
c a forward and backward substitution step using the previously-computed
c lu factors.
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
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
      call second(t1)
      do 20 i=2,nstop
      loc=i
      iord=nodplc(irswpf+i)
   10 loc=nodplc(jcpt+loc)
      if (nodplc(jcolno+loc).ge.i) go to 20
      j=nodplc(jcolno+loc)
      jord=nodplc(irswpf+j)
      value(lvn+iord)=value(lvn+iord)-
     1           value(lvn+loc)*value(lvn+jord)
      go to 10
   20 continue
c
c     back substitution
c
      i=nstop
      iord=nodplc(irswpf+i)
      jord=nodplc(icswpf+i)
      locnn=indxx(iord,jord)
   30 value(lvn+iord)=value(lvn+iord)/value(lvn+locnn)
      i=i-1
      if (i.le.1) go to 100
      iord=nodplc(irswpf+i)
      loc=i
   35 loc=nodplc(jcpt+loc)
   40 if (nodplc(jcolno+loc).ne.i) go to 35
      locnn=loc
   50 loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 30
      j=nodplc(jcolno+loc)
      jord=nodplc(irswpf+j)
      value(lvn+iord)=value(lvn+iord)-
     1           value(lvn+loc)*value(lvn+jord)
      go to 50
  100 call second(t2)
      rstats(46)=rstats(46)+t2-t1
      return
      end
