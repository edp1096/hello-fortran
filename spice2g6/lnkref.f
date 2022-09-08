      subroutine lnkref
      implicit double precision (a-h,o-z)
c
c     this routine resolves all unsatisfied name references.
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
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c  mutual inductors
c
      loc=locate(4)
  100 if (loc.eq.0) go to 200
      iref=nodplc(loc+2)
      call fndnam(value(iunsat+iref),loc-1,loc+2,3)
      iref=nodplc(loc+3)
      call fndnam(value(iunsat+iref),loc-1,loc+3,3)
      loc=nodplc(loc)
      go to 100
c
c  current-controlled current source
c
  200 loc=locate(7)
  210 if (loc.eq.0) go to 300
      nump=nodplc(loc+4)
      locp=nodplc(loc+6)
      do 220 i=1,nump
      iref=nodplc(locp+i)
      call fndnam(value(iunsat+iref),loc-1,locp+i,9)
  220 continue
      loc=nodplc(loc)
      go to 210
c
c  current-controlled voltage sources
c
  300 loc=locate(8)
  310 if (loc.eq.0) go to 400
      nump=nodplc(loc+4)
      locp=nodplc(loc+7)
      do 320 i=1,nump
      iref=nodplc(locp+i)
      call fndnam(value(iunsat+iref),loc-1,locp+i,9)
  320 continue
      loc=nodplc(loc)
      go to 310
c
c  diodes
c
  400 loc=locate(11)
  410 if (loc.eq.0) go to 500
      iref=nodplc(loc+5)
      call fndnam(value(iunsat+iref),loc-1,loc+5,21)
      loc=nodplc(loc)
      go to 410
c
c  bjts
c
  500 loc=locate(12)
  510 if (loc.eq.0) go to 600
      iref=nodplc(loc+8)
      call fndnam(value(iunsat+iref),loc-1,loc+8,22)
      loc=nodplc(loc)
      go to 510
c
c  jfets
c
  600 loc=locate(13)
  610 if (loc.eq.0) go to 700
      iref=nodplc(loc+7)
      call fndnam(value(iunsat+iref),loc-1,loc+7,23)
      loc=nodplc(loc)
      go to 610
c
c  mosfets
c
  700 loc=locate(14)
  710 if (loc.eq.0) go to 1000
      iref=nodplc(loc+8)
      call fndnam(value(iunsat+iref),loc-1,loc+8,24)
      loc=nodplc(loc)
      go to 710
c
c  finished
c
 1000 call clrmem(iunsat)
      return
      end
