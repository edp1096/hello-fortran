      subroutine subckt
      implicit double precision (a-h,o-z)
c
c     this routine drives the expansion of subcircuit calls.
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
c
c... avoid 'call by value' problems, make inodi, inodx arrays
c... in routines which receive them as parameters ]]]
      locx=locate(19)
   10 if (locx.eq.0) go to 300
      locs=nodplc(locx+3)
      asnam=value(iunsat+locs)
      call fndnam(asnam,locx-1,locx+3,20)
      if (nogo.ne.0) go to 300
      locs=nodplc(locx+3)
c
c  check for recursion
c
      isbptr=nodplc(locx-1)
   20 if (isbptr.eq.0) go to 30
      if (locs.eq.nodplc(isbptr+3)) go to 260
      isbptr=nodplc(isbptr-1)
      go to 20
c
c
   30 call sizmem(nodplc(locx+2),nxnod)
      call sizmem(nodplc(locs+2),nssnod)
      if (nxnod.ne.nssnod) go to 250
      call getm4(inodx,nssnod)
      call getm4(inodi,nssnod)
      itemp=nodplc(locs+2)
      call copy4(nodplc(itemp+1),nodplc(inodx+1),nssnod)
      itemp=nodplc(locx+2)
      call copy4(nodplc(itemp+1),nodplc(inodi+1),nxnod)
c
c  add elements of subcircuit to nominal circuit
c
      loc=nodplc(locs+3)
  100 if (loc.eq.0) go to 200
      id=nodplc(loc-1)
      if (id.eq.20) go to 110
      call find(dble(jelcnt(id)),id,loce,1)
      nodplc(loce-1)=locx
      call addelt(loce,loc,id,inodx,inodi,nxnod)
  110 loc=nodplc(loc)
      go to 100
c
c
  200 call clrmem(inodx)
      call clrmem(inodi)
      locx=nodplc(locx)
      go to 10
c
c  errors
c
  250 locv=nodplc(locx+1)
      axnam=value(locv)
      locv=nodplc(locs+1)
      asnam=value(locv)
      write (iofile,251) axnam,asnam
  251 format('0*error*:  ',a8,' has different number of nodes than ',a8/
     1)
      nogo=1
      go to 300
  260 locsv=nodplc(locs+1)
      asnam=value(locsv)
      write (iofile,261) asnam
  261 format('0*error*:  subcircuit ',a8,' is defined recursively'/)
      nogo=1
c
c  finished
c
  300 return
      end
