c spice version 2g.6  sccsid=dctran.ma 3/15/83
      subroutine dctran
      implicit double precision (a-h,o-z)
c
c
c     this routine controls the dc transfer curve, dc operating point,
c and transient analyses.  the variables mode and modedc (defined below)
c determine exactly which analysis is performed.
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
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      logical memptr
c
c
      dimension subtit(4,2)
      dimension avhdr(3),avfrm(4)
      data avhdr / 8h( (2x,a4, 8h,3x,a7,3, 5hx)//) /
      data avfrm / 8h( (1h ,a, 8h1,i3,1h), 8h,f10.4,3, 4hx)/) /
      data anode, avltg / 4hnode, 7hvoltage /
      data subtit / 8hsmall si, 8hgnal bia, 8hs soluti, 8hon      ,
     1              8hinitial , 8htransien, 8ht soluti, 8hon      /
      data lprn /1h(/
      data ablnk, aletr, alett /1h , 1hr, 1ht /
c
c      the variables *mode*, *modedc*, and *initf* are used by spice to
c keep track of the state of the analysis.  the values of these flags
c (and the corresponding meanings) are as follows:
c
c        flag    value    meaning
c        ----    -----    -------
c
c        mode      1      dc analysis (subtype defined by *modedc*)
c                  2      transient analysis
c                  3      ac analysis (small signal)
c
c        modedc    1      dc operating point
c                  2      initial operating point for transient analysis
c                  3      dc transfer curve computation
c
c        initf     1      converge with 'off' devices allowed to float
c                  2      initialize junction voltages
c                  3      converge with 'off' devices held 'off'
c                  4      store small-signal parameters away
c                  5      first timepoint in transient analysis
c                  6      prediction step
c
c note:  *modedc* is only significant if *mode* = 1.
c
c
c  initialize
c
      call second(t1)
      sfactr=1.0d0
c.. don't take any chances with lx3, set to large number
      lx3=20000000
      lx2=20000000
c.. see if lx3 and lx2 tables are needed
      nolx2=0
      nolx3=0
   20 loctim=5
c
c.. post-processing initialization
c
      if(ipostp.eq.0) go to 25
      numcur=jelcnt(9)
      numpos=nunods+numcur
      call getm8(ibuff,numpos)
      numpos=numpos*4
      if(numcur.eq.0) go to 25
      loc=locate(9)
      loccur=nodplc(loc+6)-1
c
c...  set up format
c
   25 nvprln=4+(lwidth-72)/19
      nvprln=min0(nvprln,ncnods-1)
      ipos=2
      call alfnum(nvprln,avfrm,ipos)
      ipos=2
      call alfnum(nvprln,avhdr,ipos)
c...  allocate storage
      if (mode.eq.2) go to 35
      need=4*nstop+nttbr+nxtrm
      call avlm8(navl)
      if(need.le.navl) go to 30
c...  not enough memory for dc operating point analysis
      write(iofile,26) need,navl
   26 format('0insufficient memory available for dc analysis.',/
     1' memory required ',i6,', memory available ',i6,'.')
      nogo=1
      go to 1100
   30 call getm8(lvnim1,nstop)
      call getm8(lvn,nstop+nttbr)
      call slpmem(lvn,nstop)
      call getm8(lx0,nxtrm)
      call getm8(lvntmp,nstop)
      if (modedc.ne.3) go to 45
   35 call getm8(lx1,nxtrm)
      if(nolx2.eq.0) call getm8(lx2,nxtrm)
      if (mode.ne.2) go to 40
      if(nolx3.eq.0) call getm8(lx3,nxtrm)
      call getm8(ltd,0)
   40 call getm8(loutpt,0)
   45 call crunch
   50 if (mode.eq.2) go to 500
      time=0.0d0
      ag(1)=0.0d0
      call sorupd
      if (modedc.eq.3) go to 300
c
c
c  ....  single point dc analysis
c
c
c  compute dc operating point
c
  100 if (itl6.gt.0) go to 105
      initf=2
      call iter8(itl1)
      rstats(6)=rstats(6)+iterno
      if (igoof.ne.0) go to 150
      go to 110
  105 call sorstp(itl6)
      rstats(6)=rstats(6)+iterno
      if (igoof.ne.0) go to 150
  110 if (modedc.ne.1) go to 120
      initf=4
      call diode
      call bjt
      call jfet
      call mosfet
c
c  print operating point
c
  120 if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 1000
      call title(-1,lwidth,1,subtit(1,modedc))
      write (iofile,avhdr) (anode,avltg,i=1,nvprln)
      write (iofile,avfrm) (lprn,nodplc(junode+i),value(lvnim1+i),
     1  i=2,ncnods)
      go to 1000
c
c  no convergence
c
  150 nogo=1
      write (iofile,151)
  151 format('1*error*:  no convergence in dc analysis'/'0last node vol'
     1   ,'tages:'/)
      write (iofile,avhdr) (anode,avltg,i=1,nvprln)
      write (iofile,avfrm) (lprn,nodplc(junode+i),value(lvnim1+i),
     1  i=2,ncnods)
      go to 1000
c
c  ....  dc transfer curves
c
  300 numout=jelcnt(41)+1
      if(ipostp.ne.0) call pheadr(atitle)
      itemp=itcelm(1)
      locs=nodplc(itemp+1)
      anam=value(locs)
      call move(anam,2,ablnk,1,7)
      irdctc=0
      irdct2=0
      itdctc=0
      itdct2=0
      if (anam.eq.aletr) irdctc=1
      if (anam.eq.alett) itdctc=1
      temval=value(locs+1)
      icvfl2=1
      if(itcelm(2).eq.0) go to 310
      itemp=itcelm(2)
      locs2=nodplc(itemp+1)
      anam=value(locs2)
      call move(anam,2,ablnk,1,7)
      if (anam.eq.aletr) irdct2=1
      if (anam.eq.alett) itdct2=1
      temv2=value(locs2+1)
      value(locs2+1)=tcstar(2)
      temp=dabs((tcstop(2)-tcstar(2))/tcincr(2))+0.5d0
      icvfl2=idint(temp)+1
      icvfl2=max0(icvfl2,1)
  310 delta=tcincr(1)
      do 320 i=1,7
      delold(i)=delta
  320 continue
      icvfl1=icvflg/icvfl2
      value(locs+1)=tcstar(1)
      if ((itdctc.ne.1).and.(itdct2.ne.1)) go to 325
      itemno=3
      if (itdctc.eq.1) value(itemps+itemno)=value(locs+1)
      if (itdct2.eq.1) value(itemps+itemno)=value(locs2+1)
      call tmpupd
  325 if (irdctc.eq.1) value(locs+1)=1.0d0/value(locs+1)
      if (irdct2.eq.1) value(locs2+1)=1.0d0/value(locs2+1)
      icalc=0
      ical2=0
      loctim=3
  340 initf=2
      call iter8(itl1)
      rstats(4)=rstats(4)+iterno
      call copy8(value(lx0+1),value(lx1+1),nxtrm)
      if(nolx2.eq.0) call copy8(value(lx0+1),value(lx2+1),nxtrm)
      if (igoof.ne.0) go to 450
      go to 360
  350 call getcje
      if ((maxtim-itime).le.limtim) go to 460
      initf=6
      call iter8(itl2)
      rstats(4)=rstats(4)+iterno
      if (igoof.ne.0) go to 340
c
c  store outputs
c
  360 call extmem(loutpt,numout)
      loco=loutpt+icalc*numout
      icalc=icalc+1
      ical2=ical2+1
      value(loco+1)=value(locs+1)
      if (irdctc.eq.1) value(loco+1)=1.0d0/value(loco+1)
      loc=locate(41)
  370 if (loc.eq.0) go to 400
      if (nodplc(loc+5).ne.0) go to 380
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      iseq=nodplc(loc+4)
      value(loco+iseq)=value(lvnim1+node1)-value(lvnim1+node2)
      loc=nodplc(loc)
      go to 370
  380 iptr=nodplc(loc+2)
      iptr=nodplc(iptr+6)
      iseq=nodplc(loc+4)
      value(loco+iseq)=value(lvnim1+iptr)
      loc=nodplc(loc)
      go to 370
c
c  increment source value
c
  400 if(ipostp.eq.0) go to 410
      value(ibuff+1)=value(locs+1)
      call copy8(value(lvnim1+2),value(ibuff+2),nunods-1)
      if(numcur.ne.0) call copy8(value(lvnim1+loccur+1),
     1  value(ibuff+nunods+1),numcur)
      call fwrite(value(ibuff+1),numpos)
  410 if (icalc.ge.icvflg) go to 490
      if(ical2.ge.icvfl1) go to 480
      if(nolx2.ne.0) go to 420
      call ptrmem(lx2,itemp)
      call ptrmem(lx1,lx2)
      go to 430
  420 call ptrmem(lx1,itemp)
  430 call ptrmem(lx0,lx1)
      call ptrmem(itemp,lx0)
      value(locs+1)=tcstar(1)+dble(ical2)*delta
      if (itdctc.ne.1) go to 440
      value(itemps+itemno-1)=value(itemps+itemno)
      value(itemps+itemno)=value(locs+1)
      call tmpupd
  440 if (irdctc.eq.1) value(locs+1)=1.0d0/value(locs+1)
      go to 350
c
c  no convergence
c
  450 itemp=itcelm(1)
      loce=nodplc(itemp+1)
      write (iofile,451) value(loce),value(locs+1)
  451 format('1*error*:  no convergence in dc transfer curves at ',a8,
     1   ' = ',1pd10.3/'0last node voltages:'/)
      write (iofile,avhdr) (anode,avltg,i=1,nvprln)
      write (iofile,avfrm) (lprn,nodplc(junode+i),value(lvnim1+i),
     1  i=2,ncnods)
      go to 470
  460 write (iofile,461)
  461 format('0*error*:  cpu time limit exceeded ... analysis stopped'/)
      go to 470
  462 write(iofile,463)
  463 format('0*error*:   temperature sweep should be the second sweep
     1source, change the order and re-execute'/)
  470 nogo=1
      go to 490
c... reset first sweep variable ... step second
  480 ical2=0
      value(locs+1)=tcstar(1)
      if (irdctc.eq.1) value(locs+1)=1.0d0/value(locs+1)
      if (itdctc.eq.1) go to 462
      value(locs2+1)=value(locs2+1)+tcincr(2)
      if (irdct2.eq.1) value(locs2+1)=1.0d0/value(locs2+1)
      if (itdct2.ne.1) go to 340
      value(itemps+itemno-1)=value(itemps+itemno)
      value(itemps+itemno)=value(locs2+1)
      call tmpupd
      go to 340
c
c  finished with dc transfer curves
c
  490 value(locs+1)=temval
      if(itcelm(2).ne.0) value(locs2+1)=temv2
      if ((itdctc.eq.0).and.(itdct2.eq.0)) go to 1000
      value(itemps+itemno-1)=value(itemps+itemno)
      if (itdctc.eq.1) value(itemps+itemno)=temval
      if (itdct2.eq.1) value(itemps+itemno)=temv2
      write (iofile,492)
  492 format (/,'0*****0 return to original temperature 0*****0',/)
      call tmpupd
      itemno=1
      call relmem(itemps,2)
      if(ipostp.eq.0) go to 1000
      call fwrite(value(ibuff+1),numpos)
      go to 1000
c
c  ....  transient analysis
c
  500 numout=jelcnt(42)+1
      if(ipostp.ne.0) call pheadr(atitle)
c...  limit delmax if no energy-storage elements
      numese=jelcnt(2)+jelcnt(3)+jelcnt(11)+jelcnt(12)+jelcnt(13)
     1   +jelcnt(14)
      if (numese.eq.0) delmax=dmin1(delmax,tstep)
      initf=5
      iord=1
      loctim=9
      icalc=0
      numtp=0
      numrtp=0
      numnit=0
      time=0.0d0
      ibkflg=1
      delbkp=delmax
      nbkpt=1
      delta=delmax
      do 510 i=1,7
      delold(i)=delta
  510 continue
      delnew=delta
      delmin=1.0d-9*delmax
      go to 650
c
c  increment time, update sources, and solve next timepoint
c
  600 time=time+delta
      call sorupd
      if (nogo.ne.0) go to 950
      call getcje
      if ((maxtim-itime).le.limtim) go to 920
      if ((itl5.ne.0).and.(numnit.ge.itl5)) go to 905
      call comcof
      if (initf.ne.5) initf=6
      itrlim=itl4
      if ((numtp.eq.0).and.(nosolv.ne.0)) itrlim=itl1
      call iter8(itrlim)
      numnit=numnit+iterno
      numtp=numtp+1
      if (numtp.ne.1) go to 605
      if(nolx2.eq.0) call copy8(value(lx1+1),value(lx2+1),nxtrm)
      if(nolx3.eq.0) call copy8(value(lx1+1),value(lx3+1),nxtrm)
c.. note that time-point is cut when itrlim exceeded regardless
c.. of which time-step contol is specified thru 'lvltim'.
  605 if (igoof.eq.0) go to 610
      jord=iord
      iord=1
      if (jord.ge.5) call clrmem(lx7)
      if (jord.ge.4) call clrmem(lx6)
      if (jord.ge.3) call clrmem(lx5)
      if ((jord.ge.2).and.(method.ne.1)) call clrmem(lx4)
      igoof=0
      time=time-delta
      delta=delta/8.0d0
      go to 620
  610 delnew=delta
      if (numtp.eq.1) go to 630
      call trunc(delnew)
      if (delnew.ge.(0.9d0*delta)) go to 630
      time=time-delta
      delta=delnew
  620 numrtp=numrtp+1
      ibkflg=0
      delold(1)=delta
      if (delta.ge.delmin) go to 600
      time=time+delta
      go to 900
c
c  determine order of integration method
c
c...  skip if trapezoidal algorithm used
  630 if ((method.eq.1).and.(iord.eq.2)) go to 650
      if (numtp.eq.1) go to 650
      ordrat=1.05d0
      if (iord.gt.1) go to 635
      iord=2
      call trunc(delnew)
      iord=1
      if ((delnew/delta).le.ordrat) go to 650
      if (maxord.le.1) go to 650
      iord=2
      if (method.eq.1) go to 650
      call getm8(lx4,nxtrm)
      go to 650
  635 if (iord.lt.maxord) go to 640
      iord=iord-1
      call trunc(delnew)
      iord=iord+1
      if ((delnew/delta).le.ordrat) go to 650
      go to 642
  640 iord=iord-1
      call trunc(delnew)
      iord=iord+1
      if ((delnew/delta).le.ordrat) go to 645
  642 iord=iord-1
      if (iord.eq.1) call clrmem(lx4)
      if (iord.eq.2) call clrmem(lx5)
      if (iord.eq.3) call clrmem(lx6)
      if (iord.eq.4) call clrmem(lx7)
      go to 650
  645 iord=iord+1
      call trunc(delnew)
      iord=iord-1
      if ((delnew/delta).le.ordrat) go to 650
      iord=iord+1
      if (iord.eq.2) call getm8(lx4,nxtrm)
      if (iord.eq.3) call getm8(lx5,nxtrm)
      if (iord.eq.4) call getm8(lx6,nxtrm)
      if (iord.eq.5) call getm8(lx7,nxtrm)
c
c  store outputs
c
  650 if ((time+delta).le.tstart) go to 685
      if ((numtp.eq.0).and.(nosolv.ne.0)) go to 685
      call extmem(loutpt,numout)
      loco=loutpt+icalc*numout
      icalc=icalc+1
      value(loco+1)=time
      loc=locate(42)
  670 if (loc.eq.0) go to 682
      if (nodplc(loc+5).ne.0) go to 680
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      iseq=nodplc(loc+4)
      value(loco+iseq)=value(lvnim1+node1)-value(lvnim1+node2)
      loc=nodplc(loc)
      go to 670
  680 iptr=nodplc(loc+2)
      iptr=nodplc(iptr+6)
      iseq=nodplc(loc+4)
      value(loco+iseq)=value(lvnim1+iptr)
      loc=nodplc(loc)
      go to 670
  682 if(ipostp.eq.0) go to 684
      value(ibuff+1)=time
      call copy8(value(lvnim1+2),value(ibuff+2),nunods-1)
      if(numcur.ne.0) call copy8(value(lvnim1+loccur+1),
     1  value(ibuff+nunods+1),numcur)
      call fwrite(value(ibuff+1),numpos)
  684 continue
c
c  update transmission line delay table
c
  685 if (jelcnt(17).eq.0) go to 694
      call sizmem(ltd,ltdsiz)
      numtd=ltdsiz/ntlin
      if (numtd.le.3) go to 689
      baktim=time-tdmax
      if (baktim.lt.0.0d0) go to 689
      lcntr=0
      ltemp=ltd
      do 686 i=1,numtd
      if (value(ltemp+1).ge.baktim) go to 687
      ltemp=ltemp+ntlin
      lcntr=lcntr+1
  686 continue
      go to 689
  687 if (lcntr.le.2) go to 689
      lcntr=lcntr-2
      nwords=lcntr*ntlin
      ltemp=ltemp-ntlin-ntlin
      call copy8(value(ltemp+1),value(ltd+1),ltdsiz-nwords)
      call relmem(ltd,nwords)
      call sizmem(ltd,ltdsiz)
  689 call extmem(ltd,ntlin)
      ltdptr=ltd+ltdsiz
      value(ltdptr+1)=time
      loc=locate(17)
  690 if (loc.eq.0) go to 693
      locv=nodplc(loc+1)
      z0=value(locv+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      ibr1=nodplc(loc+8)
      ibr2=nodplc(loc+9)
      lspot=nodplc(loc+30)+ltdptr
      if ((initf.eq.5).and.(nosolv.ne.0)) go to 691
      value(lspot)=value(lvnim1+node3)-value(lvnim1+node4)
     1   +value(lvnim1+ibr2)*z0
      value(lspot+1)=value(lvnim1+node1)-value(lvnim1+node2)
     1   +value(lvnim1+ibr1)*z0
      go to 692
  691 value(lspot)=value(locv+7)+value(locv+8)*z0
      value(lspot+1)=value(locv+5)+value(locv+6)*z0
  692 loc=nodplc(loc)
      go to 690
c
c  add two *fake* backpoints to ltd for interpolation near time=0.0d0
c
  693 if (numtd.ne.0) go to 694
      call extmem(ltd,ntlin+ntlin)
      call copy8(value(ltd+1),value(ltd+ntlin+1),ntlin)
      call copy8(value(ltd+1),value(ltd+2*ntlin+1),ntlin)
      value(ltd+2*ntlin+1)=time
      value(ltd+ntlin+1)=time-delta
      value(ltd+1)=time-delta-delta
c
c  rotate state vector storage
c
c.. time-point accepted
  694 call copy8(delold(1),delold(2),6)
      delta=delnew
      delold(1)=delta
      go to (710,706,702,698,696,696), iord
  696 call ptrmem(lx7,itemp)
      call ptrmem(lx6,lx7)
      go to 700
  698 call ptrmem(lx6,itemp)
  700 call ptrmem(lx5,lx6)
      go to 704
  702 call ptrmem(lx5,itemp)
  704 call ptrmem(lx4,lx5)
      go to 708
  706 if (method.eq.1) go to 710
      call ptrmem(lx4,itemp)
  708 call ptrmem(lx3,lx4)
      go to 713
  710 if(nolx3.eq.0) go to 712
      if(nolx2.eq.0) go to 711
      call ptrmem(lx1,itemp)
      go to 714
  711 call ptrmem(lx2,itemp)
      call ptrmem(lx1,lx2)
      go to 714
  712 call ptrmem(lx3,itemp)
  713 call ptrmem(lx2,lx3)
      call ptrmem(lx1,lx2)
  714 call ptrmem(lx0,lx1)
      call ptrmem(itemp,lx0)
c
c  check breakpoints
c
  750 if (ibkflg.eq.0) go to 760
c.. just accepted analysis at breakpoint
      jord=iord
      iord=1
      if (jord.ge.5) call clrmem(lx7)
      if (jord.ge.4) call clrmem(lx6)
      if (jord.ge.3) call clrmem(lx5)
      if ((jord.ge.2).and.(method.ne.1)) call clrmem(lx4)
      ibkflg=0
      nbkpt=nbkpt+1
      if (nbkpt.gt.numbkp) go to 950
      temp=dmin1(delbkp,value(lsbkpt+nbkpt)-time)
      delta=dmin1(delta,0.1d0*temp,delmax)
      if (numtp.eq.0) delta=delta/10.0d0
      delold(1)=delta
      go to 600
  760 del1=value(lsbkpt+nbkpt)-time
      if ((1.01d0*delta).le.del1) go to 600
      ibkflg=1
      delbkp=delta
      delta=del1
      delold(1)=delta
      go to 600
c
c  transient analysis failed
c
  900 write (iofile,901)
  901 format('1*error*:  internal timestep too small in transient analys
     1is'/)
      go to 910
  905 write (iofile,906) itl5
  906 format('1*error*:  transient analysis iterations exceed limit of '
     1,i5,/'0this limit may be overridden using the itl5 parameter on th
     2e .option card')
  910 write (iofile,911) time,delta,numnit
  911 format(1h0,10x,'time = ',1pd12.5,';  delta = ',d12.5,';  numnit =
     1',i6/)
      write (iofile,916)
  916 format(1h0/'0last node voltages:'/)
      write (iofile,avhdr) (anode,avltg,i=1,nvprln)
      write (iofile,avfrm) (lprn,nodplc(junode+i),value(lvnim1+i),
     1  i=2,ncnods)
      go to 930
  920 write (iofile,921) time
  921 format('0*error*:  cpu time limit exceeded in transient analysis '
     1   ,'at time = ',1pd13.6/)
  930 nogo=1
c
c  finished with transient analysis
c
  950 rstats(10)=rstats(10)+numnit
      rstats(30)=rstats(30)+numtp
      rstats(31)=rstats(31)+numrtp
      rstats(32)=rstats(32)+numnit
      if(ipostp.eq.0) go to 1000
      if (ipostp.ne.0) call clsraw
c
c  return unneeded memory
c
 1000 if (mode.eq.2) go to 1010
      if (modedc.ne.3) go to 1100
 1010 call clrmem(lvnim1)
      call clrmem(lx0)
      call clrmem(lvn)
      call clrmem(lx1)
      if (memptr(macins)) call clrmem(macins)
      if(nolx2.eq.0) call clrmem(lx2)
      call clrmem(lvntmp)
      if ((mode.eq.1).and.(modedc.eq.3)) go to 1020
      if(nolx3.eq.0) call clrmem(lx3)
      if (mode.eq.1) go to 1020
      call clrmem(ltd)
      if (iord.eq.1) go to 1020
      if (method.eq.1) go to 1020
      call clrmem(lx4)
      if (iord.eq.2) go to 1020
      call clrmem(lx5)
      if (iord.eq.3) go to 1020
      call clrmem(lx6)
      if (iord.eq.4) go to 1020
      call clrmem(lx7)
 1020 call extmem(loutpt,2*numout)
 1100 if(ipostp.ne.0) call clrmem(ibuff)
      call second(t2)
      rstats(loctim)=rstats(loctim)+t2-t1
      return
      end
