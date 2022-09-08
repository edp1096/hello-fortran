      subroutine fndnam(anam,jsbptr,ispot,id)
      implicit double precision (a-h,o-z)
c
c     this routine searches for an element with id 'id' by tracing back
c up the subcircuit definition list.  if the element is not found, the
c nominal element list is searched.
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
      integer xxor
c
c
      isbptr=nodplc(jsbptr)
   10 if (isbptr.eq.0) go to 50
      isub=nodplc(isbptr+3)
      loc=nodplc(isub+3)
   20 if (loc.eq.0) go to 40
      if (id.ne.nodplc(loc-1)) go to 30
      locv=nodplc(loc+1)
      if (xxor(anam,value(locv)).ne.0) go to 30
      if (id.ne.20) go to 50
      go to 65
   30 loc=nodplc(loc)
      go to 20
   40 isbptr=nodplc(isbptr-1)
      go to 10
c
   50 loc=locate(id)
   60 if (loc.eq.0) go to 90
      if (nodplc(loc-1).ne.isbptr) go to 70
      locv=nodplc(loc+1)
      if (xxor(anam,value(locv)).ne.0) go to 70
   65 nodplc(ispot)=loc
      go to 100
   70 loc=nodplc(loc)
      go to 60
   90 write (iofile,91) anam
   91 format('0*error*:  unable to find ',a8/)
      nogo=1
  100 return
      end
