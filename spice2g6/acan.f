c spice version 2g.6  sccsid=acan.ma 3/15/83
      subroutine acan
      implicit double precision (a-h,o-z)
c
c     this routine drives the small-signal analyses.
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
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      complex cendor
      equivalence (value(1),nodplc(1),cvalue(1))
      call second(t1)
c.. post-processor initialization
      if(ipostp.eq.0) go to 1
      numcur=jelcnt(9)
      numpos=nunods+numcur
      call getm16(ibuff,numpos)
      numpos=numpos*4
      if(numcur.eq.0) go to 1
      loc=locate(9)
      loccur=nodplc(loc+6)-1
c
c  allocate storage
c
    1 call getm8(ndiag,2*nstop)
      call getm8(lvn,nstop+nttbr)
      call getm8(imvn,nstop+nttbr)
      call getm16(lcvn,nstop)
      if (idist.ne.0) call dinit
      nandd=0
      if (inoise.eq.0) go to 10
      if (idist.eq.0) go to 10
      nandd=1
      call getm16(lvntmp,nstop)
   10 call getm16(loutpt,0)
      call crunch
      numout=jelcnt(43)+jelcnt(44)+jelcnt(45)+1
      lynl=lvn
      imynl=imvn
      lcvntp=lvntmp
      icalc=0
      if (ipostp.ne.0) call pheadr(atitle)
      freq=fstart
c
c  load y matrix and c vector, solve for v vector
c
  100 call getcje
      if ((maxtim-itime).le.limtim) go to 900
      omega=twopi*freq
      call acload
  110 call acdcmp
      call acsol
      if (igoof.eq.0) go to 200
      write (iofile,121) igoof,freq
  121 format('0warning:  underflow ',i4,' time(s) in ac analysis at freq
     1 = ',1pd9.3,' hz')
      igoof=0
c
c  store outputs
c
  200 call extmem(loutpt,numout)
      loco=loutpt+icalc*numout
      icalc=icalc+1
      cvalue(loco+1)=cmplx(sngl(freq),sngl(omega))
      loc=locate(43)
  310 if (loc.eq.0) go to 350
      if (nodplc(loc+5).ne.0) go to 320
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      iseq=nodplc(loc+4)
      cvalue(loco+iseq)=cvalue(lcvn+node1)-cvalue(lcvn+node2)
      loc=nodplc(loc)
      go to 310
  320 iptr=nodplc(loc+2)
      iptr=nodplc(iptr+6)
      iseq=nodplc(loc+4)
      cvalue(loco+iseq)=cvalue(lcvn+iptr)
      loc=nodplc(loc)
      go to 310
  350 if(ipostp.eq.0) go to 400
      cvalue(ibuff+1)=cmplx(sngl(freq),0.0e0)
      call copy16(cvalue(lcvn+2),cvalue(ibuff+2),nunods-1)
      if(numcur.ne.0) call copy16(cvalue(lcvn+loccur+1),
     1  cvalue(ibuff+nunods+1),numcur)
      call dblsgl(cvalue(ibuff+1),numpos)
      call fwrite(cvalue(ibuff+1),numpos)
c
c  noise and distortion analyses
c
  400 if (nandd.eq.0) go to 410
      call copy16(cvalue(lcvn+1),cvalue(lcvntp+1),nstop)
  410 if (inoise.ne.0) call noise(loco)
      if (nandd.eq.0) go to 420
      call copy16(cvalue(lcvntp+1),cvalue(lcvn+1),nstop)
  420 if (idist.ne.0) call disto(loco)
c
c  increment frequency
c
      if (icalc.ge.jacflg) go to 1000
      if (idfreq.ge.3) go to 510
      freq=freq*fincr
      go to 100
  510 freq=freq+fincr
      go to 100
c
c  finished
c
  900 write (iofile,901)
  901 format('0*error*:  cpu time limit exceeded ... analysis stopped'/)
      nogo=1
 1000 if(ipostp.eq.0) go to 1010
      if (ipostp.ne.0) call clsraw
      if(ipostp.ne.0) call clrmem(ibuff)
 1010 call clrmem(lvnim1)
      call clrmem(lx0)
      call clrmem(lvn)
      call clrmem(imvn)
      call clrmem(lcvn)
      call clrmem(ndiag)
      if (idist.eq.0) go to 1020
      call clrmem(ld0)
      call clrmem(ld1)
 1020 if (nandd.eq.0) go to 1040
      call clrmem(lvntmp)
 1040 call second(t2)
      rstats(7)=rstats(7)+t2-t1
      rstats(8)=rstats(8)+icalc
      return
      end
