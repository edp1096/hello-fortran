      subroutine acdcmp
      implicit double precision (a-h,o-z)
c
c     this routine performs an lu factorization of the circuit equation
c coefficient matrix.
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
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
      n=1
   10 n=n+1
      nxti=n
      nxtj=n
c
c     calculate contribution from (nxti,nxtj)
c
      if (n.ge.nstop) return
      n1=nodplc(irswpf+nxti)
      n2=nodplc(icswpf+nxtj)
      locnn=indxx(n1,n2)
      gdiag=dabs(value(lynl+locnn))+dabs(value(imynl+locnn))
      if (gdiag.ge.pivtol) go to 20
      value(lynl+locnn)=pivtol
      value(imynl+locnn)=0.0d0
      write(iofile,11) n
   11 format(1h0,' underflow occured at step n= ',i5)
c
c     down col j
c
   20 locr=nodplc(irpt+locnn)
   25 if (locr.eq.0) go to 10
      i=nodplc(irowno+locr)
      call cdiv(value(lynl+locr),value(imynl+locr),value(lynl+locnn),
     1     value(imynl+locnn),value(lynl+locr),value(imynl+locr))
      locc=nodplc(jcpt+locnn)
c
c     for each element look up row nxti
c
   30 if (locc.eq.0) go to 70
      j=nodplc(jcolno+locc)
c
c     locate element (i,j)
c
   35 if (j.lt.i) go to 45
      locij=locc
   40 locij=nodplc(irpt+locij)
      if (nodplc(irowno+locij).eq.i) go to 55
      go to 40
   45 locij=locr
   50 locij=nodplc(jcpt+locij)
      if (nodplc(jcolno+locij).eq.j) go to 55
      go to 50
   55 call cmult(value(lynl+locc),value(imynl+locc),
     1     value(lynl+locr),value(imynl+locr),xreal,ximag)
      value(lynl+locij)=value(lynl+locij)-xreal
      value(imynl+locij)=value(imynl+locij)-ximag
      locc=nodplc(jcpt+locc)
      go to 30
   70 locr=nodplc(irpt+locr)
      go to 25
      end
