      subroutine dcdcmp
      implicit double precision (a-h,o-z)
c
c     this routine swaps rows and columns in the coefficient matrix accor-
c ding to the numerical pivoting and minimum fillin terms.it then performs
c an in-place lu factorization of the coefficient matrix.
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
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      call second(t1)
      if (ipiv.le.0) go to 12
      if (idebug(11).le.0) go to  3
      call dmpmat(6hdcdcmp)
      idebug(11)=idebug(11)-1
    3 do 10 i=2,nstop
      no=0
      loc=nodplc(jcpt+i)
    5 if (loc.eq.0) go to 7
      no=no+1
      loc=nodplc(jcpt+loc)
      go to 5
    7 nodplc(numoff+i)=no
      no=0
      loc=nodplc(irpt+i)
    8 if (loc.eq.0) go to 9
      no=no+1
      loc=nodplc(irpt+loc)
      go to 8
    9 nodplc(nmoffc+i)=no
   10 continue
   12 n=1
c
c     find next pivot
c
   15 n=n+1
      if (ipiv.gt.0) go to 20
      if (idebug(13).le.0) go to 17
      call dmpmat(6hdcdcm2)
      idebug(13)=idebug(13)-1
   17 if (idebug(14).le.0) go to 18
      if (mode.ne.2) go to 18
      call dmpmat(6hdcdcm3)
      idebug(14)=idebug(14)-1
   18 if (idebug(15).le.0.or.icalc.le.10) go to 19
      call dmpmat(6hdcdcm4)
      idebug(15)=idebug(15)-1
   19 nxti=n
      nxtj=n
      go to 120
c
c     search the coresponding column for max entry
c
   20 vmax=0.0d0
      loci=n
   25 loci=nodplc(irpt+loci)
      if (loci.eq.0) go to 50
      i=nodplc(irowno+loci)
      if (i.lt.n) go to 25
   30 if (dabs(value(lvn+loci)).le.vmax) go to 25
      vmax=dabs(value(lvn+loci))
      go to 25
   50 if (vmax.gt.pivtol) go to 60
      write(iofile,51) n,vmax
   51 format('0*error*:  maximum entry in this column at step ',i4,' (',
     1   1pd12.6,') is less than pivtol')
      igoof=1
      return
   60 epsrel=dmax1(pivrel*vmax,pivtol)
      if (n.ge.nstop) go to 200
      if (ipiv.le.0) go to 120
c
c     pivoting on the diagonal
c
      minop=100000
      nxti=0
      do 70 i=n,nstop
      i1=nodplc(irswpf+i)
      j1=nodplc(icswpf+i)
      ispot=indxx(i1,j1)
      if (ispot.eq.1) go to 70
      if (dabs(value(lvn+ispot)).lt.epsrel) go to 70
      nop=(nodplc(numoff+i)-1)*(nodplc(nmoffc+i)-1)
      if (nop.ge.minop) go to 70
      minop=nop
      nxti=i
      nxtj=i
      if (minop.le.0) go to 95
   70 continue
      if (nxti.le.0) go to 75
      if (nxti-n) 120,120,100
c
c     pivoting on the entire matrix
c
   75 do 90 i=n,nstop
      loc=i
   80 loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 90
      j=nodplc(jcolno+loc)
      if (j.lt.n) go to 80
      if (dabs(value(lvn+loc)).lt.epsrel) go to 80
      nop=(nodplc(numoff+i)-1)*(nodplc(nmoffc+j)-1)
      if (nop.ge.minop) go to 80
      minop=nop
      nxti=i
      nxtj=j
      if (minop.le.0) go to 95
   90 continue
      if (nxti.gt.0) go to 95
      write (iofile,92)
   92 format('0*abort*:  pivot not in dcdcmp')
      igoof=1
      go to 200
   95 if (nxti.eq.n.and.nxtj.eq.n) go to 120
      if (nxti.eq.n) go to 105
c
c     a pivot has been found
c
  100 load=nodplc(irswpf+nxti)
      lr=nodplc(irswpf+n)
      nodplc(irswpf+nxti)=lr
      nodplc(irswpr+lr)=nxti
      nodplc(irswpf+n)=load
      nodplc(irswpr+load)=n
      noff=nodplc(numoff+nxti)
      nodplc(numoff+nxti)=nodplc(numoff+n)
      nodplc(numoff+n)=noff
      if (nxtj.eq.n) go to 110
  105 load=nodplc(icswpf+nxtj)
      lc=nodplc(icswpf+n)
      nodplc(icswpf+nxtj)=lc
      nodplc(icswpr+lc)=nxtj
      nodplc(icswpf+n)=load
      nodplc(icswpr+load)=n
      noff=nodplc(nmoffc+nxtj)
      nodplc(nmoffc+nxtj)=nodplc(nmoffc+n)
      nodplc(nmoffc+n)=noff
  110 call swapij(nxti,n,nxtj,n)
      nxti=n
      nxtj=n
c
c     calculate contribution from nxti, nxtj and find fill-ins
c
  120 if (n.ge.nstop) go to 200
      n1=nodplc(irswpf+nxti)
      n2=nodplc(icswpf+nxtj)
      locnn=indxx(n1,n2)
      if (ipiv.le.0 .and. dabs(value(lvn+locnn)).lt.pivtol) go to 220
c
c     down col j
c
      locr=nodplc(irpt+locnn)
  125 if (locr.eq.0) go to 180
      i=nodplc(irowno+locr)
      value(lvn+locr)=value(lvn+locr)/value(lvn+locnn)
      locc=nodplc(jcpt+locnn)
c
c     for each column element look up row nxti
c
  130 if (locc.eq.0) go to 170
      j=nodplc(jcolno+locc)
c
c     check for fill-in (i,j)
c
      if (ipiv.le.0) go to 135
      call sizmem(jcpt,isize1)
      call reserv(i,j)
      call sizmem(jcpt,isize2)
      if (isize1.eq.isize2) go to 135
      call extmem(lvn,1)
      nttbr=nttbr+1
      value(lvn+nstop+nttbr)=0.0d0
c
c     locate element (i,j)
c
  135 if (j.lt.i) go to 145
      locij=locc
  140 locij=nodplc(irpt+locij)
      if (nodplc(irowno+locij).eq.i) go to 155
      go to 140
  145 locij=locr
  150 locij=nodplc(jcpt+locij)
      if (nodplc(jcolno+locij).eq.j) go to 155
      go to 150
  155 value(lvn+locij)=value(lvn+locij)-
     1                  value(lvn+locc)*value(lvn+locr)
  160 locc=nodplc(jcpt+locc)
      go to 130
  170 locr=nodplc(irpt+locr)
      if (ipiv.le.0) go to 125
      nodplc(numoff+i)=nodplc(numoff+i)-1
      go to 125
c
c     reduce nmoffc for each element in col nxti
c
  180 if (ipiv.le.0) go to 15
      locc=nodplc(jcpt+locnn)
  185 if (locc.eq.0) go to 15
      j=nodplc(jcolno+locc)
      nodplc(nmoffc+j)=nodplc(nmoffc+j)-1
      locc=nodplc(jcpt+locc)
      go to 185
c
c     done
c
  200 if (ipiv.eq.0) go to 210
      if (idebug(17).le.0) go to 202
      call dmpmat(6hdcdcm5)
      idebug(17)=idebug(17)-1
  202 call matloc
      rstats(49)=rstats(49)+1.0d0
      ipiv=0
      if (lvlcod.eq.2) lvlcod=3
      ifill=nttbr-nttar
      perspa=100.0d0*(1.0d0-dble(nttbr)/dble(nstop*nstop))
c
c  calculation of operation count (operation := `*' or `/'):
c
c     noffr := off-diagonal elements in row, not including diagonal,
c                counting only those elements in the remainder matrix
c     noffc := off-diagonal elements in column, not including diagonal,
c                counting only those elements in the remainder matrix
c
c     then we have
c
c        lu decomposition     requires sigma(2,nstop-1) {noffc + noffc*noffr}
c        forward substitution          sigma(2,nstop-1) {noffc + 1}   +   1
c        backward substitution         sigma(2,nstop-1) {noffr}
c
c     which sums to
c
c               sigma(2,nstop-1) {noffc + noffc*noffr + (noffc+1) + noffr} + 1
c         or
c               sigma(2,nstop-1) {noffc*(noffr+2) + noffr + 1}   +   1
c
      iops=1
      nstop1=nstop-1
      do 205 i=2,nstop1
      noffr=nodplc(numoff+i)-1
      noffc=nodplc(nmoffc+i)-1
      iops=iops+noffr+noffc*(noffr+2)+1
  205 continue
      rstats(20)=nstop
      rstats(21)=nttar
      rstats(22)=nttbr
      rstats(23)=ifill
      rstats(24)=0.0d0
      rstats(25)=nttbr
      rstats(26)=iops
      rstats(27)=perspa
      go to 215
  210 if (idebug(18).le.0) go to 212
      call dmpmat(6hdcdcm6)
      idebug(18)=idebug(18)-1
  212 if (idebug(19).le.0.or.icalc.le.10) go to 215
      call dmpmat(6hdcdcm7)
      idebug(19)=idebug(19)-1
  215 call second(t2)
      rstats(50)=rstats(50)+t2-t1
       return
  220 ipiv=1
      write(iofile,221) n,nxti,nxtj,iterno,time
  221 format(' pivot change on fly: n= ',i5,' nxti= ',i5,' nxtj= ',
     1       i5,' iterno= ',i5,' time= ',1pd12.5)
      rstats(49)=rstats(49)+1.0d0
      go to 20
      end
