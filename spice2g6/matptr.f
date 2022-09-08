      subroutine matptr
      implicit double precision (a-h,o-z)
c
c     this routine (by calls to the routine reserve) establishes the
c nonzero-element structure of the circuit equation coefficient matrix.
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
c  allocate and initialize storage
c
      call getm4(isr,nstop+1)
      numvs=jelcnt(3)+jelcnt(6)+jelcnt(8)+jelcnt(9)+2*jelcnt(17)
      call getm4(iseq,numvs)
      call getm4(iseq1,numvs)
      call getm4(neqn,numvs)
      call getm4(nodevs,numnod)
      call getm4(ndiag,nstop)
      call getm4(nmoffc,nstop)
      call getm4(numoff,nstop)
      call getm4(irpt,nstop)
      call getm4(jcpt,nstop)
      call getm4(irowno,nstop)
      call getm4(jcolno,nstop)
      call slpmem(irpt,nstop)
      call slpmem(jcpt,nstop)
      call slpmem(irowno,nstop)
      call slpmem(jcolno,nstop)
      call crunch
c
      call zero4(nodplc(irpt+1),nstop)
      call zero4(nodplc(jcpt+1),nstop)
      call zero4(nodplc(irowno+1),nstop)
      call zero4(nodplc(jcolno+1),nstop)
      call zero4(nodplc(iseq1+1),numvs)
      call zero4(nodplc(nodevs+1),numnod)
      call zero4(nodplc(ndiag+1),nstop)
      call zero4(nodplc(nmoffc+1),nstop)
      call zero4(nodplc(numoff+1),nstop)
c
      numvs=0
      nxtrm=0
      ndist=0
      ntlin=1
      ibr=numnod
c
c  resistors
c
      loc=locate(1)
  110 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 120
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      call reserv(node1,node1)
      call reserv(node1,node2)
      call reserv(node2,node1)
      call reserv(node2,node2)
      loc=nodplc(loc)
      go to 110
c
c  capacitors
c
  120 loc=locate(2)
  130 if ((loc.eq.0).or.(nodplc(loc+12).ne.0)) go to 400
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      call reserv(node1,node2)
      call reserv(node2,node1)
      ntemp=nodplc(ndiag+node1)
      call reserv(node1,node1)
      nodplc(ndiag+node1)=ntemp
      ntemp=nodplc(ndiag+node2)
      call reserv(node2,node2)
      nodplc(ndiag+node2)=ntemp
      nodplc(loc+8)=nxtrm+1
      nxtrm=nxtrm+2
      loc=nodplc(loc)
      go to 130
c
c  inductors
c
  400 loc=locate(3)
  430 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 440
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ibr=ibr+1
      nodplc(loc+5)=ibr
      call reserv(node1,ibr)
      call reserv(node2,ibr)
      call reserv(ibr,node1)
      call reserv(ibr,node2)
      ntemp=nodplc(ndiag+ibr)
      call reserv(ibr,ibr)
      nodplc(ndiag+ibr)=ntemp
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(neqn+numvs)=ibr
      nodplc(nodevs+node1)=nodplc(nodevs+node1)+1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)+1
      nodplc(loc+11)=nxtrm+1
      nxtrm=nxtrm+2
      loc=nodplc(loc)
      go to 430
c
c  mutual inductors
c
  440 loc=locate(4)
  450 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 460
      nl1=nodplc(loc+2)
      nl2=nodplc(loc+3)
      nl1=nodplc(nl1+5)
      nl2=nodplc(nl2+5)
      call reserv(nl1,nl2)
      call reserv(nl2,nl1)
      loc=nodplc(loc)
      go to 450
c
c  nonlinear voltage-controlled current sources
c
  460 loc=locate(5)
  462 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 464
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      ndim2=ndim+ndim
      locn=nodplc(loc+6)
      do 463 i=1,ndim2
      node=nodplc(locn+i)
      call reserv(node1,node)
      call reserv(node2,node)
  463 continue
      nodplc(loc+12)=nxtrm+1
      nxtrm=nxtrm+1+ndim2
      loc=nodplc(loc)
      go to 462
c
c  nonlinear voltage controlled voltage sources
c
  464 loc=locate(6)
  466 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 468
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ibr=ibr+1
      nodplc(loc+6)=ibr
      call reserv(node1,ibr)
      call reserv(node2,ibr)
      call reserv(ibr,node1)
      call reserv(ibr,node2)
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(neqn+numvs)=ibr
      nodplc(nodevs+node1)=nodplc(nodevs+node1)+1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)+1
      ndim=nodplc(loc+4)
      ndim2=ndim+ndim
      locn=nodplc(loc+7)
      do 467 i=1,ndim2
      node=nodplc(locn+i)
      call reserv(ibr,node)
  467 continue
      nodplc(loc+13)=nxtrm+1
      nxtrm=nxtrm+2+ndim2
      loc=nodplc(loc)
      go to 466
c
c  voltage sources
c
  468 loc=locate(9)
  470 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 472
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ibr=ibr+1
      nodplc(loc+6)=ibr
      call reserv(node1,ibr)
      call reserv(node2,ibr)
      call reserv(ibr,node1)
      call reserv(ibr,node2)
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(neqn+numvs)=ibr
      nodplc(nodevs+node1)=nodplc(nodevs+node1)+1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)+1
      loc=nodplc(loc)
      go to 470
c
c  nonlinear current controlled current sources
c
  472 loc=locate(7)
  474 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 476
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      locvs=nodplc(loc+6)
      do 475 i=1,ndim
      locvst=nodplc(locvs+i)
      kbr=nodplc(locvst+6)
      call reserv(node1,kbr)
      call reserv(node2,kbr)
  475 continue
      nodplc(loc+12)=nxtrm+1
      nxtrm=nxtrm+1+ndim+ndim
      loc=nodplc(loc)
      go to 474
c
c  nonlinear current controlled voltage sources
c
  476 loc=locate(8)
  478 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 500
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ibr=ibr+1
      nodplc(loc+6)=ibr
      call reserv(node1,ibr)
      call reserv(node2,ibr)
      call reserv(ibr,node1)
      call reserv(ibr,node2)
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(neqn+numvs)=ibr
      nodplc(nodevs+node1)=nodplc(nodevs+node1)+1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)+1
      ndim=nodplc(loc+4)
      locvs=nodplc(loc+7)
      do 479 i=1,ndim
      locvst=nodplc(locvs+i)
      kbr=nodplc(locvst+6)
      call reserv(ibr,kbr)
  479 continue
      nodplc(loc+13)=nxtrm+1
      nxtrm=nxtrm+2+ndim+ndim
      loc=nodplc(loc)
      go to 478
c
c  diodes
c
  500 loc=locate(11)
  510 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 520
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      call reserv(node1,node1)
      call reserv(node2,node2)
      call reserv(node3,node3)
      call reserv(node1,node3)
      call reserv(node2,node3)
      call reserv(node3,node1)
      call reserv(node3,node2)
      nodplc(loc+11)=nxtrm+1
      nxtrm=nxtrm+5
      nodplc(loc+12)=ndist+1
      ndist=ndist+7
      loc=nodplc(loc)
      go to 510
c
c  transistors
c
  520 loc=locate(12)
  530 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 540
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      node7=nodplc(loc+30)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      cdis=value(locm+32)
      call reserv(node1,node1)
      call reserv(node2,node2)
      call reserv(node3,node3)
      call reserv(node4,node4)
      call reserv(node5,node5)
      call reserv(node6,node6)
      call reserv(node1,node4)
      call reserv(node2,node5)
      call reserv(node3,node6)
      call reserv(node4,node5)
      call reserv(node4,node6)
      call reserv(node5,node6)
      call reserv(node4,node1)
      call reserv(node5,node2)
      call reserv(node6,node3)
      call reserv(node5,node4)
      call reserv(node6,node4)
      call reserv(node6,node5)
      call reserv(node7,node7)
      call reserv(node4,node7)
      call reserv(node7,node4)
      if (cdis.lt.1.0d0) call reserv(node2,node4)
      if (cdis.lt.1.0d0) call reserv(node4,node2)
      nodplc(loc+22)=nxtrm+1
      nxtrm=nxtrm+19
      nodplc(loc+23)=ndist+1
      ndist=ndist+21
      loc=nodplc(loc)
      go to 530
c
c  jfets
c
  540 loc=locate(13)
  550 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 560
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      call reserv(node1,node1)
      call reserv(node2,node2)
      call reserv(node3,node3)
      call reserv(node4,node4)
      call reserv(node5,node5)
      call reserv(node1,node4)
      call reserv(node2,node4)
      call reserv(node2,node5)
      call reserv(node3,node5)
      call reserv(node4,node5)
      call reserv(node4,node1)
      call reserv(node4,node2)
      call reserv(node5,node2)
      call reserv(node5,node3)
      call reserv(node5,node4)
      nodplc(loc+19)=nxtrm+1
      nxtrm=nxtrm+13
      loc=nodplc(loc)
      go to 550
c
c  mosfets
c
  560 loc=locate(14)
  570 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 600
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      call reserv(node1,node1)
      call reserv(node2,node2)
      call reserv(node3,node3)
      call reserv(node4,node4)
      call reserv(node5,node5)
      call reserv(node6,node6)
      call reserv(node1,node5)
      call reserv(node2,node4)
      call reserv(node2,node5)
      call reserv(node2,node6)
      call reserv(node3,node6)
      call reserv(node4,node5)
      call reserv(node4,node6)
      call reserv(node5,node6)
      call reserv(node5,node1)
      call reserv(node4,node2)
      call reserv(node5,node2)
      call reserv(node6,node2)
      call reserv(node6,node3)
      call reserv(node5,node4)
      call reserv(node6,node4)
      call reserv(node6,node5)
      nodplc(loc+26)=nxtrm+1
      nxtrm=nxtrm+28
      loc=nodplc(loc)
      go to 570
c
c  transmission lines
c
  600 loc=locate(17)
  610 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 1000
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      ni1=nodplc(loc+6)
      ni2=nodplc(loc+7)
      ibr1=ibr+1
      ibr2=ibr+2
      ibr=ibr+2
      nodplc(loc+8)=ibr1
      nodplc(loc+9)=ibr2
      call reserv(node1,node1)
      call reserv(node1,ni1)
      call reserv(node2,ibr1)
      call reserv(node3,node3)
      call reserv(node4,ibr2)
      call reserv(ni1,node1)
      call reserv(ni1,ni1)
      call reserv(ni1,ibr1)
      call reserv(ni2,ni2)
      call reserv(ni2,ibr2)
      call reserv(ibr1,node2)
      call reserv(ibr1,node3)
      call reserv(ibr1,node4)
      call reserv(ibr1,ni1)
      call reserv(ibr1,ibr2)
      call reserv(ibr2,node1)
      call reserv(ibr2,node2)
      call reserv(ibr2,node4)
      call reserv(ibr2,ni2)
      call reserv(ibr2,ibr1)
      call reserv(node3,ni2)
      call reserv(ni2,node3)
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(iseq1+numvs)=1
      nodplc(neqn+numvs)=ibr1
      nodplc(nodevs+ni1)=nodplc(nodevs+ni1)+1
      nodplc(nodevs+node2)=nodplc(nodevs+node2)+1
      numvs=numvs+1
      nodplc(iseq+numvs)=loc
      nodplc(iseq1+numvs)=2
      nodplc(neqn+numvs)=ibr2
      nodplc(nodevs+ni2)=nodplc(nodevs+ni2)+1
      nodplc(nodevs+node4)=nodplc(nodevs+node4)+1
      nodplc(loc+30)=ntlin+1
      ntlin=ntlin+2
      loc=nodplc(loc)
      go to 610
c
c  finished
c
 1000 return
      end
