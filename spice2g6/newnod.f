      subroutine newnod(nodold,nodnew,inodx,inodi,nnodi)
      implicit double precision (a-h,o-z)
c
c     this routine makes a new node number for an element which is about
c to be added to the circuit as a result of a subcircuit call.
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
c... inodx, inodi are arrays (see subckt)
      dimension inodx(1),inodi(1)
c
      if (nodold.ne.0) go to 5
      nodnew=1
      go to 20
    5 do 10 i=1,nnodi
      jnodx=inodx(1)
      if (nodold.ne.nodplc(jnodx+i)) go to 10
      jnodi=inodi(1)
      nodnew=nodplc(jnodi+i)
      go to 20
   10 continue
c
      call extmem(inodx(1),1)
      call extmem(inodi(1),1)
      call extmem(junode,1)
      nnodi=nnodi+1
      ncnods=ncnods+1
      jnodx=inodx(1)
      nodplc(jnodx+nnodi)=nodold
      jnodi=inodi(1)
      nodplc(jnodi+nnodi)=ncnods
      nodplc(junode+ncnods)=nodplc(junode+ncnods-1)+1
      nodnew=ncnods
   20 return
      end
