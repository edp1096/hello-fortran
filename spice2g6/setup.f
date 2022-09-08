c spice version 2g.6  sccsid=setup.ma 3/15/83
      subroutine setup
      implicit double precision (a-h,o-z)
c
c     this routine drives the sparse matrix setup used by spice.
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      logical memptr
c
      call second(t1)
      nstop=numnod+jelcnt(3)+jelcnt(6)+jelcnt(8)+jelcnt(9)+2*jelcnt(17)
c
c     clear old tables
c
      if (memptr(irpt)) call clrmem(irpt)
      if (memptr(jcpt)) call clrmem(jcpt)
      if (memptr(irowno)) call clrmem(irowno)
      if (memptr(jcolno)) call clrmem(jcolno)
c
c  reserve matrix locations for each element
c
      call matptr
      if (nogo.ne.0) go to 1000
c
c  reorder matrix pointers
c
      nttbr=0
      do 120 i=2,nstop
      loc=i
  110 if (nodplc(jcpt+loc).eq.0) go to 120
      loc=nodplc(jcpt+loc)
      nttbr=nttbr+1
      go to 110
  120 continue
c...  add ground
      nttar=nttbr
      call reordr
      if (nogo.ne.0) go to 1000
c
c  store matrix locations
c
      call matloc
c
c  .nodeset
c
  200 call sizmem(nsnod,nic)
      if(nic.eq.0) go to 220
      call getm4(nsmat,nic)
      do 210 i=1,nic
      node=nodplc(nsnod+i)
      nodplc(nsmat+i)=indxx(node,node)
  210 continue
c
c  transient initial conditions
c
  220 call sizmem(icnod,nic)
      if(nic.eq.0) go to 300
      call getm4(icmat,nic)
      do 230 i=1,nic
      node=nodplc(icnod+i)
      nodplc(icmat+i)=indxx(node,node)
  230 continue
c
  300 call clrmem(iseq)
      call clrmem(iseq1)
      call clrmem(neqn)
      call clrmem(nodevs)
      call clrmem(ndiag)
c
c  finished
c
 1000 call second(t2)
      rstats(2)=rstats(2)+t2-t1
      return
      end
