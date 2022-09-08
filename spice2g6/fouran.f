      subroutine fouran
      implicit double precision (a-h,o-z)
c
c     this routine determines the fourier coefficients of a transient
c analysis waveform.
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=outinf 3/15/83
      common /outinf/ xincr,string(15),xstart,yvar(8),itab(8),itype(8),
     1   ilogy(8),npoint,numout,kntr,numdgt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension sinco(9),cosco(9)
      dimension fortit(4)
      data fortit / 8hfourier , 8hanalysis, 8h        , 8h         /
      data ablnk / 1h  /
c
c
      forprd=1.0d0/forfre
      xstart=tstop-forprd
      kntr=1
cc    xn=101.0d0
      xincr=forprd/npoint
cc    npoint=xn
      call getm8(locx,npoint)
      call getm8(locy,npoint)
      do 105 nknt=1,nfour
      itab(1)=nodplc(ifour+nknt)
      kfrout=itab(1)
      call ntrpl8(locx,locy,numpnt)
      dcco=0.0d0
      call zero8(sinco,9)
      call zero8(cosco,9)
      loct=locy+1
      ipnt=0
   10 yvr=value(loct+ipnt)
      dcco=dcco+yvr
      forfac=dble(ipnt)*twopi/npoint
      arg=0.0d0
      do 20 k=1,9
      arg=arg+forfac
      sinco(k)=sinco(k)+yvr*dsin(arg)
      cosco(k)=cosco(k)+yvr*dcos(arg)
   20 continue
      ipnt=ipnt+1
      if (ipnt.ne.npoint) go to 10
      dcco=dcco/npoint
      forfac=2.0d0/npoint
      do 30 k=1,9
      sinco(k)=sinco(k)*forfac
      cosco(k)=cosco(k)*forfac
   30 continue
      call title(0,72,1,fortit)
      ipos=1
      call outnam(kfrout,1,string,ipos)
      call move(string,ipos,ablnk,1,7)
      jstop=(ipos+6)/8
      write (iofile,61) (string(j),j=1,jstop)
   61 format(' fourier components of transient response ',5a8///)
      write (iofile,71) dcco
   71 format('0dc component =',1pd12.3/,
     1   '0harmonic   frequency    fourier    normalized    phase     no
     2rmalized'/,
     3   '    no         (hz)     component    component    (deg)    pha
     4se (deg)'//)
      iknt=1
      freq1=forfre
      xnharm=1.0d0
      call magphs(cmplx(sngl(sinco(1)),sngl(cosco(1))),xnorm,pnorm)
      phasen=0.0d0
      write (iofile,81) iknt,freq1,xnorm,xnharm,pnorm,phasen
   81 format(i6,1pd15.3,d12.3,0pf13.6,f10.3,f12.3/)
      thd=0.0d0
      do 90 iknt=2,9
      freq1=dble(iknt)*forfre
      call magphs(cmplx(sngl(sinco(iknt)),sngl(cosco(iknt))),
     1   harm,phase)
      xnharm=harm/xnorm
      phasen=phase-pnorm
      thd=thd+xnharm*xnharm
      write (iofile,81) iknt,freq1,harm,xnharm,phase,phasen
   90 continue
      thd=100.0d0*dsqrt(thd)
      write (iofile,101) thd
  101 format (//5x,'total harmonic distortion =  ',f12.6,'  percent')
  105 continue
      call clrmem(locx)
      call clrmem(locy)
  110 return
      end
