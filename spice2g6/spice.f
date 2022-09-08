c spice version 2g.6  sccsid=root.ma 3/15/83
      program spice
c...  cray notes:
c.. change data filnam to 5lspice
c.. change call overlay .. remove last zero
c.. change ??????program to subroutine in all but this overlay
c.. delete overlay spice,0,0 card below
c..
      implicit double precision (a-h,o-z)
c
c
c spice version 2g.6  sccsid=root.ma 3/15/83
c
c     spice is an electronic circuit simulation program that was deve-
c loped by the integrated circuits group of the electronics research
c laboratory and the department of electrical engineering and computer
c sciences at the university of california, berkeley, california.  the
c program spice is available free of charge to any interested party.
c the sale, resale, or use of this program for profit without the
c express written consent of the department of electrical engineering
c and computer sciences, university of california, berkeley, california,
c is forbidden.
c
c
c     the present version is based on the spice2 program versions 2e.3
c and 2f.1 developed at the university of california berkeley and the
c hewlett-packard spice version 2.7.
c     this version is designed to be transportable on most computers
c with an ansi fortran compiler and enough memory space for code and
c data. the memory manager uses the function 'locf' to find the
c address of a pointer; this function must be provided.
c
c
c implementation notes:
c
c     subroutines mclock and mdate return the time (as hh:mm:ss) and
c the date (as dd mmm yy), respectively.  subroutine getcje returns in
c common block /cje/ various attributes of the current job environment.
c spice expects getcje to set /cje/ variables maxtim, itime, and icost.
c maxtim is the maximum cpu time in seconds, itime is the elapsed cpu
c time in seconds, and icost is the job cost in cents.
c subroutine memory is used to change the number of memory words
c allocated to spice.  if the amount of memory allocated to a jobstep
c is fixed, subroutine memory need not be changed.
c     subroutine second(t) returns the time in seconds and is used
c for timing purposes. it must be provided where not available.
c     ifamwa (set in a data statement below) should be set to the
c address of the first available word of memory (following overlays, if
c any).  the proper value should be easily obtainable from any load map.
c ifamwa is used only on computers where the program (spice) can change
c the allocated memory dynamically at run time according to circuit size.
c (see also comments under subroutine setmem).
c     all berkeley spice2.f release versions do not implement the ifamwa
c feature due to its dependence on operating system.
c     with the exception of most flags, all data in spice are stored in
c the form of managed tables allocated in the /blank/ array value().
c array value() can be redimensioned in the main program according to
c memory availability at each user site. it should be noted again that
c the program dynamically manages its data within the bounds of array
c value().
c     the vax release versions assume the virtual memory feature and
c dimension value() to 200,000 double precision words.
c     the cdc and ibm versions dimension value() to 20000 real or
c real*8 words, respectively.
c     spice is particularly well-suited to being run using a one-level
c overlay structure beginning with routines spice (the overlay root),
c readin, errchk, setup, dctran, dcop, acan, and ovtpvt.  the order of
c the routines in this listing corresponds to that structure.  note
c that if cdc explicit overlaying is to be used, an overlay directive
c card must be inserted before the first line of each of the just-named
c routines.
c
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
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
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
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=outinf 3/15/83
      common /outinf/ xincr,string(15),xstart,yvar(8),itab(8),itype(8),
     1   ilogy(8),npoint,numout,kntr,numdgt
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
c
      dimension acctit(4)
      dimension remain(4)
      data ablnk /1h  /
      data acctit / 8hjob stat, 8histics s, 8hummary  , 8h         /
      data ahdr1,ahdr2,ahdr3 / 8h  spice ,8h2g.6    ,8h3/15/83     /
c
c
      ipostp=0
c     check if a raw data file should be written
c     open file if so required
      if  (iargc().gt.0) ipostp=iopraw()
      maxmem=2*200000
      maxtim=1e8
      icost=0
      iofile=6
c
c  initialization
c
      aprog(1)=ahdr1
      aprog(2)=ahdr2
      aprog(3)=ahdr3
      achar=ablnk
      keof=0
      call xtime(atime)
      call xdate(adate)
      boltz=1.3806226d-23
      charge=1.6021918d-19
      ctok=273.15d0
      eps0=8.854214871d-12
      epssil=11.7d0*eps0
      epsox=3.9d0*eps0
      twopi=8.0d0*datan2(1.0d0,1.0d0)
      rad=360.0d0/twopi
      xlog2=dlog(2.0d0)
      xlog10=dlog(10.0d0)
      root2=dsqrt(2.0d0)
      nodata=1
c
c  begin job
c
   10 if (keof.eq.1) go to 1000
      call getcje
      call second(REAL(time1))
      icost1=icost
      igoof=0
      mode=0
      nogo=0
      call setmem(nodplc(1),maxmem)
      if (nogo.ne.0) go to 1000
      call zero8(rstats,50)
c
c  read remainder of data deck and check for input errors
c
      call readin
      if (nogo.ne.0) go to 300
      if (keof.eq.1) go to 1000
      nodata=0
   50 call errchk
      if (nogo.ne.0) go to 300
      call setup
      if (nogo.ne.0) go to 300
c
c  change parameters and re-analisis
c
      if (numcyc.eq.0) go to 90
   70 if (numcyc.ge.numalt) go to 310
      numcyc=numcyc+1
      call alter
      if (nogo.ne.0) go to 305
c
c  cycle through temperatures
c
   90 itemno=1
      if (numtem.eq.1) go to 110
  100 if (itemno.eq.numtem) go to 70
      itemno=itemno+1
      call tmpupd
c
c  dc transfer curves
c
  110 if (icvflg.eq.0) go to 150
c...  see routine *dctran* for explanation of *mode*, etc.
      mode=1
      modedc=3
      call dctran
      call ovtpvt
      if (nogo.ne.0) go to 300
c
c  small signal operating point
c
  150 if (kssop.gt.0) go to 170
      if (jacflg.ne.0) go to 170
      if ((icvflg+jtrflg).gt.0) go to 250
  170 mode=1
      modedc=1
      call dctran
      if (nogo.ne.0) go to 300
      call dcop
      if (nogo.ne.0) go to 300
c
c  ac small signal analysis
c
  200 if (jacflg.eq.0) go to 250
      mode=3
      call acan
      call ovtpvt
      if (nogo.ne.0) go to 300
c
c  transient analysis
c
  250 if (jtrflg.eq.0) go to 100
      mode=1
      modedc=2
      call dctran
      if (nogo.ne.0) go to 300
      call dcop
      if (nogo.ne.0) go to 300
      mode=2
      call dctran
      call ovtpvt
      if (nogo.ne.0) go to 300
      go to 100
c
c  job concluded
c
  300 write (iofile,301)
  301 format(1h0,9x,'***** job aborted')
      nodata=0
      go to 320
  305 write (iofile,306)
  306 format (1h0,9x,'***** this parameter change is illegal')
  310 write (iofile,311)
  311 format(1h0,/,9x,'job concluded')
c
c  job accounting
c
  320 continue
      numel=0
      do 360 i=1,18
  360 numel=numel+jelcnt(i)
      numtem=max0(numtem-1,1)
      idist=min0(idist,1)
      if (iprnta.eq.0) go to 800
      call title(-1,lwidth,1,acctit)
      write (iofile,361) nunods,ncnods,numnod,numel,(jelcnt(i),i=11,14)
  361 format('   nunods ncnods numnod numel  diodes  bjts  jfets  mfets'
     1   //,i9,2i7,i6,i8,i6,2i7)
      write (iofile,371) numtem,icvflg,jtrflg,jacflg,inoise,idist,nogo
  371 format(/'0  numtem icvflg jtrflg jacflg inoise  idist   nogo'/,
     1   2h0 ,7i7)
      write (iofile,381) rstats(20),rstats(21),rstats(22),rstats(23),
     1   rstats(26),rstats(27)
  381 format(/'0  nstop   nttbr   nttar   ifill    iops    perspa'//,
     1   1x,5f8.0,f9.3)
      write (iofile,391) rstats(30),rstats(31),rstats(32),maxmem,maxuse,
     1   cpyknt
  391 format(/'0  numttp  numrtp  numnit  maxmem  memuse  copyknt',//,
     1   2x,3f8.0,2x,i6,2x,i6,2x,f8.0)
      write (iofile,401) (rstats(i),i=1,6),rstats(50),rstats(49),
     1   rstats(46),(rstats(i),i=7,11)
  401 format(/,
     1   1h0,9x,'readin  ',12x,f10.2/,
     2   1h0,9x,'setup   ',12x,f10.2/,
     3   1h0,9x,'trcurv  ',12x,f10.2,10x,f6.0/,
     4   1h0,9x,'dcan    ',12x,f10.2,10x,f6.0/,
     5   1h0,9x,'dcdcmp  ',12x,f10.3,10x,f6.0/,
     6   1h0,9x,'dcsol   ',12x,f10.3/,
     7   1h0,9x,'acan    ',12x,f10.2,10x,f6.0/,
     8   1h0,9x,'tranan  ',12x,f10.2,10x,f6.0/,
     9   1h0,9x,'output  ',12x,f10.2)
      write (iofile,402) rstats(45),rstats(48),rstats(47),rstats(44),
     1   rstats(43)
  402 format(
     1   1h0,9x,'load    ',12x,f10.3/,
     2   1h0,9x,'codgen  ',12x,f10.3,10x,f6.0/,
     3   1h0,9x,'codexc  ',12x,f10.3/,
     4   1h0,9x,'macins  ',12x,f10.3)
  800 call getcje
      call second(REAL(time2))
      et=time2-time1
      tcost=dble(icost-icost1)/100.0d0
      if (iprnta.eq.0) go to 810
      ohead=et-(rstats(1)+rstats(2)+rstats(3)+rstats(5)+rstats(7)
     1   +rstats(9)+rstats(11))
      write (iofile,801) ohead
  801 format(1h0,9x,'overhead',12x,f10.2)
  810 write (iofile,811) et
  811 format(1h0,9x,'total job time      ',f10.2)
      rstats(33)=cpyknt
      rstats(34)=et
      rstats(35)=tcost
      rstats(36)=ohead
  900 if ((maxtim-itime).ge.limtim) go to 10
      write (iofile,901)
  901 format('1warning:  further analysis stopped due to cpu time limit'
     1/)
 1000 if(nodata.ne.0) write(iofile,1001)
 1001 format(/1x,'input deck (file) contains no data.')
      stop
      end
