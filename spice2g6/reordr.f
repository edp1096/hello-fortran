      subroutine reordr
      implicit double precision (a-h,o-z)
c
c     this routine swaps rows in the coefficient matrix to eliminate
c singularity problems which can be recognized by examining the circuit
c topology.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c  allocate and initialize storage
c
      call getm4(irswpf,nstop)
      call getm4(irswpr,nstop)
      call getm4(icswpf,nstop)
      call getm4(icswpr,nstop)
c
      do 10 i=1,nstop
      nodplc(irswpf+i)=i
   10 continue
      call copy4(nodplc(irswpf+1),nodplc(irswpr+1),nstop)
      call copy4(nodplc(irswpf+1),nodplc(icswpf+1),nstop)
      call copy4(nodplc(irswpf+1),nodplc(icswpr+1),nstop)
c
c  swap current equations into admittance part of equation matrix
c
      nextv=1
c
c  find suitable voltage source
c
  100 if (nextv.gt.numvs) go to 600
      ix=0
      do 130 i=nextv,numvs
      loc=nodplc(iseq+i)
      node=nodplc(loc+2)
      nflag=nodplc(iseq1+i)
      if (nflag.eq.1) node=nodplc(loc+6)
      if (nflag.eq.2) node=nodplc(loc+7)
      if (node.eq.1) go to 110
      if (nodplc(nodevs+node).ge.2) go to 110
      if (nodplc(ndiag+node).eq.0) go to 145
      ix=i
      locx=loc
      nodex=node
  110 node=nodplc(loc+3)
      if (nflag.eq.2) node=nodplc(loc+5)
      if (node.eq.1) go to 130
      if (nodplc(nodevs+node).ge.2) go to 130
  120 if (nodplc(ndiag+node).eq.0) go to 145
      ix=i
      locx=loc
      nodex=node
  130 continue
      if (ix.eq.0) go to 590
      i=ix
      loc=locx
      node=nodex
c
c  resequence voltage sources
c
  145 nodplc(iseq+i)=nodplc(iseq+nextv)
      nodplc(iseq+nextv)=loc
      ltemp=nodplc(iseq1+i)
      nodplc(iseq1+i)=nodplc(iseq1+nextv)
      nodplc(iseq1+nextv)=ltemp
      ibr=nodplc(neqn+i)
      nodplc(neqn+i)=nodplc(neqn+nextv)
      nodplc(neqn+nextv)=ibr
      node1=nodplc(loc+2)
      if (ltemp.eq.1) node1=nodplc(loc+6)
      if (ltemp.eq.2) node1=nodplc(loc+7)
      node2=nodplc(loc+3)
      if (ltemp.eq.1) node2=nodplc(loc+3)
      if (ltemp.eq.2) node2=nodplc(loc+5)
      nodplc(nodevs+node1)=nodplc(nodevs+node1)-1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)-1
c
c  set row swap indicators
c
      l=nodplc(irswpf+ibr)
      j=nodplc(irswpr+node)
      nodplc(irswpf+j)=l
      nodplc(irswpr+l)=j
      nodplc(irswpf+ibr)=node
      nodplc(irswpr+node)=ibr
      call swapij(ibr,j,1,1)
      nextv=nextv+1
      go to 100
c
c
c  error - voltage-source/inductor/transmission-line loop detected ...
c
  590 nogo=1
      write (iofile,591)
c...  loop should have been detected in topchk
  591 format('0*abort*:  spice internal error in reordr'/)
c
c  finished
c
  600 return
      end
