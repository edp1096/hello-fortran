      subroutine topchk
      implicit double precision (a-h,o-z)
c
c     this routine constructs the element node table.  it also checks
c for voltage source/inductor loops, current source/capacitor cutsets,
c and that every node has a dc (conductive) path to ground
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
      integer change
c
c
      dimension atable(12),aide(20),nnods(20)
      dimension idlist(4),idlis2(4)
      dimension toptit(4)
      data toptit / 8helement , 8hnode tab, 8hle      , 8h         /
      data idlist / 3, 6, 8, 9 /
      data idlis2 /14,14,14,11 /
      data aide / 1hr,0.0d0,1hl,2*0.0d0,1he,0.0d0,1hh,1hv,0.0d0,1hd,
     1   1hq,1hj,1hm,0.0d0,0.0d0,1ht,0.0d0,0.0d0,0.0d0 /
      data nnods / 2,2,2,0,2,2,2,2,2,2,2,4,3,4,4,4,4,0,0,0 /
      data ablnk /1h /
c
c  allocate storage
c
      call getm4(iorder,ncnods)
      call getm4(iur,ncnods+1)
c
c  construct node table
c
      kntlim=lwidth/11
 1300 call getm4(itable,0)
      call getm4(itabid,0)
      istop=ncnods+1
      do 1310 i=1,istop
 1310 nodplc(iur+i)=1
      do 1370 id=1,18
      if (nnods(id).eq.0) go to 1370
      loc=locate(id)
 1320 if (loc.eq.0) go to 1370
      nloc=loc+1
      jstop=nnods(id)
 1330 do 1360 j=1,jstop
      node=nodplc(nloc+j)
      ispot=nodplc(iur+node+1)
      k=nodplc(iur+ncnods+1)
      call extmem(itable,1)
      call extmem(itabid,1)
      if (k.le.ispot) go to 1340
      call copy4(nodplc(itable+ispot),nodplc(itable+ispot+1),k-ispot)
      call copy4(nodplc(itabid+ispot),nodplc(itabid+ispot+1),k-ispot)
 1340 nodplc(itable+ispot)=loc
      nodplc(itabid+ispot)=id
c...  treat the substrate node of a mosfet as if it were a transmission
c...  line node, i.e. let it dangle if desired
      if(id.eq.14.and.j.eq.4) nodplc(itabid+ispot)=17
      k=node
      kstop=ncnods+1
 1350 k=k+1
      if (k.gt.kstop) go to 1360
      nodplc(iur+k)=nodplc(iur+k)+1
      go to 1350
 1360 continue
      loc=nodplc(loc)
      go to 1320
 1370 continue
c
c  check that every node has a dc path to ground
c
      call zero4(nodplc(iorder+1),ncnods)
      nodplc(iorder+1)=1
 1420 iflag=0
      do 1470 i=2,ncnods
      if (nodplc(iorder+i).eq.1) go to 1470
      jstart=nodplc(iur+i)
      jstop=nodplc(iur+i+1)-1
      if (jstart.gt.jstop) go to 1470
      do 1450 j=jstart,jstop
      loc=nodplc(itable+j)
      id=nodplc(itabid+j)
      if (aide(id).eq.0.0d0) go to 1450
      if (id.eq.17) go to 1445
      kstop=loc+nnods(id)-1
      do 1440 k=loc,kstop
      node=nodplc(k+2)
      if (nodplc(iorder+node).eq.1) go to 1460
 1440 continue
      go to 1450
 1445 if (nodplc(loc+2).eq.i) node=nodplc(loc+3)
      if (nodplc(loc+3).eq.i) node=nodplc(loc+2)
      if (nodplc(loc+4).eq.i) node=nodplc(loc+5)
      if (nodplc(loc+5).eq.i) node=nodplc(loc+4)
      if (nodplc(iorder+node).eq.1) go to 1460
 1450 continue
      go to 1470
 1460 nodplc(iorder+i)=1
      iflag=1
 1470 continue
      if (iflag.eq.1) go to 1420
c
c  print node table and topology error messages
c
      if (iprntn.eq.0) go to 1510
      call title(0,lwidth,1,toptit)
 1510 do 1590 i=1,ncnods
      jstart=nodplc(iur+i)
      jstop=nodplc(iur+i+1)-1
      if (iprntn.eq.0) go to 1550
      if (jstart.le.jstop) go to 1520
      write (iofile,1511) nodplc(junode+i)
 1511 format(1h0,i7)
      go to 1550
 1520 kntr=0
      jflag=1
      do 1540 j=jstart,jstop
      loc=nodplc(itable+j)
      locv=nodplc(loc+1)
      kntr=kntr+1
      atable(kntr)=value(locv)
      if (kntr.lt.kntlim) go to 1540
      if (jflag.eq.0) go to 1525
      jflag=0
      write (iofile,1521) nodplc(junode+i),(atable(k),k=1,kntr)
 1521 format(1h0,i7,3x,12(1x,a8))
      go to 1530
 1525 write (iofile,1526) (atable(k),k=1,kntr)
 1526 format(11x,12(1x,a8))
 1530 kntr=0
 1540 continue
      if (kntr.eq.0) go to 1550
      if (jflag.eq.0) go to 1545
      write (iofile,1521) nodplc(junode+i),(atable(k),k=1,kntr)
      go to 1550
 1545 write (iofile,1526) (atable(k),k=1,kntr)
 1550 if (jstart-jstop) 1560,1552,1556
c
c  allow node with only one connection iff element is a t-line
c
 1552 if (nodplc(itabid+jstart).eq.17) go to 1560
 1556 nogo=1
      write (iofile,1557) nodplc(junode+i)
 1557 format('0*error*:  less than 2 connections at node ',i6/)
      go to 1590
 1560 if (nodplc(iorder+i).eq.1) go to 1590
      nogo=1
      write (iofile,1561) nodplc(junode+i)
 1561 format('0*error*:  no dc path to ground from node ',i6/)
 1590 continue
c
c  check for inductor/voltage source loops
c
      do 1700 i=1,ncnods
      call zero4(nodplc(iorder+1),ncnods)
      nodplc(iorder+i)=-1
 1605 change=0
      do 1690 idcntr=1,4
      id=idlist(idcntr)
      loc=locate(id)
 1610 if ((loc.eq.0).or.(nodplc(loc+idlis2(idcntr)).ne.0)) go to 1690
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      if (nodplc(iorder+node1).eq.loc.or.
     1   nodplc(iorder+node2).eq.loc) go to 1680
      if (nodplc(iorder+node1)) 1620,1640,1630
 1620 nodplc(iorder+node1)=loc
      change=1
 1630 node=node2
      go to 1670
 1640 if (nodplc(iorder+node2)) 1650,1680,1660
 1650 nodplc(iorder+node2)=loc
      change=1
 1660 node=node1
 1670 if (nodplc(iorder+node).ne.0) go to 1710
      nodplc(iorder+node)=loc
      change=1
 1680 loc=nodplc(loc)
      go to 1610
 1690 continue
      if (change.eq.1) go to 1605
 1700 continue
      go to 1900
c ... loop found
 1710 locv=nodplc(loc+1)
      write (iofile,1711) value(locv)
 1711 format('0*error*:  inductor/voltage source loop found, containing
     1',a8/)
      nogo=1
c
c
 1900 call clrmem(iorder)
      call clrmem(iur)
      call clrmem(itable)
      call clrmem(itabid)
 2000 return
      end
