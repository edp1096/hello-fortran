      subroutine iter8(itlim)
      implicit double precision (a-h,o-z)
c
c     this routine drives the newton-raphson iteration technique used to
c solve the set of nonlinear circuit equations.
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      igoof=0
      iterno=0
      ndrflo=0
      noncon=0
      ipass=0
c
c  construct linear equations and check convergence
c
   10 ivmflg=0
      call load
   15 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 300
      iterno=iterno+1
      go to (20,30,40,60,50,60),initf
   20 if(mode.ne.1) go to 22
      call sizmem(nsnod,nic)
      if (nic.eq.0) go to 22
      if (ipass.ne.0) noncon=ipass
      ipass=0
   22 if (noncon.eq.0) go to 300
      go to 100
   30 initf=3
      if(lvlcod.eq.3) lvlcod=2
      ipiv=1
   40 if (noncon.eq.0) initf=1
      ipass=1
      go to 100
   50 if (iterno.gt.1) go to 60
      ipiv=1
      if (lvlcod.eq.3) lvlcod=2
   60 initf=1
c
c  solve equations for next iteration
c
  100 if (iterno.ge.itlim) go to 200
  102 call dcdcmp
      if (igoof.ne.0) go to 400
      if (lvlcod.eq.1) go to 105
  105 call dcsol
      go to 120
  120 if (igoof.eq.0) go to 130
      igoof=0
      if (lvlcod.ne.1) lvlcod=2
      ipiv=1
      call load
      go to 102
  130 value(lvn+1)=0.0d0
      do 135 i=1,nstop
      j=nodplc(icswpr+i)
      k=nodplc(irswpf+j)
      value(lvntmp+k)=value(lvnim1+i)
  135 continue
      call copy8(value(lvntmp+1),value(lvnim1+1),nstop)
      ntemp=noncon
      noncon=0
      if (ntemp.gt.0) go to 150
      if (iterno.eq.1) go to 150
      do 140 i=2,numnod
      vold=value(lvnim1+i)
      vnew=value(lvn+i)
      tol=reltol*dmax1(dabs(vold),dabs(vnew))+vntol
      if (dabs(vold-vnew).le.tol) go to 140
      noncon=noncon+1
  140 continue
  150 do 160 i=1,nstop
      j=nodplc(icswpr+i)
      k=nodplc(irswpf+j)
      value(lvnim1+i)=value(lvn+k)
  160 continue
c     write(iofile,151) (value(lvn+k),k=1,nstop)
c 151 format(' solution: '/1p12d10.3)
      go to 10
c
c  no convergence
c
  200 igoof=1
  300 if (ndrflo.eq.0) go to 400
      write (iofile,301) ndrflo
  301 format('0warning:  underflow occurred ',i4,' time(s)')
c
c  finished
c
  400 return
      end
