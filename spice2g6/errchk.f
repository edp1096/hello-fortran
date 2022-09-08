c spice version 2g.6  sccsid=errchk.ma 3/15/83
      subroutine errchk
      implicit double precision (a-h,o-z)
c
c
c     this routine drives the pre-processing and general error-checking
c of input performed by spice.
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
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension titlop(4)
      dimension nnods(50),aname(2)
      data aname / 4htrap, 4hgear /
      data titlop / 8hoption s, 8hummary  , 8h        , 8h         /
      data ndefin / 2h.u /
      data nnods / 2, 2, 2, 0, 2, 2, 2, 2, 2, 2,
     1             2, 4, 3, 4, 0, 0, 4, 0, 1, 0,
     2             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4             2, 2, 2, 0, 0, 0, 0, 0, 0, 0 /
      data aelmt,amodel,aoutpt /7helement,5hmodel,6houtput/
      data alsdc,alstr,alsac / 2hdc, 4htran, 2hac /
c
c
      call second(t1)
      do 60 id=1,50
      loc=locate(id)
   10 if (loc.eq.0) go to 60
      if (nodplc(loc+2).ne.ndefin) go to 50
      nogo=1
      locv=nodplc(loc+1)
      if (id.ge.21) go to 20
      anam=aelmt
      go to 40
   20 if (id.ge.31) go to 30
      anam=amodel
      go to 40
   30 anam=aoutpt
   40 write (iofile,41) anam,value(locv)
   41 format('0*error*:  ',2a8,' has been referenced but not defined'/)
   50 loc=nodplc(loc)
      go to 10
   60 continue
      if (nogo.ne.0) go to 2000
c
c  construct ordered list of user specified nodes
c
      call getm4(junode,1)
      nodplc(junode+1)=0
      nunods=1
      do 180 id=1,50
      if (nnods(id).eq.0) go to 180
      loc=locate(id)
  110 if (loc.eq.0) go to 180
      if (id.le.4) go to 120
      if (id.le.8) go to 150
      if (id.eq.19) go to 165
      if (id.le.40) go to 120
      if (id.le.43) go to 170
  120 jstop=loc+nnods(id)-1
      do 130 j=loc,jstop
      call putnod(nodplc(j+2))
  130 continue
      go to 170
  150 call putnod(nodplc(loc+2))
      call putnod(nodplc(loc+3))
      if (id.ge.7) go to 170
      locp=nodplc(loc+id+1)
      nssnod=2*nodplc(loc+4)
  155 do 160 j=1,nssnod
      call putnod(nodplc(locp+j))
  160 continue
      go to 170
  165 locp=nodplc(loc+2)
      call sizmem(nodplc(loc+2),nssnod)
      go to 155
  170 loc=nodplc(loc)
      go to 110
  180 continue
      if (nogo.ne.0) go to 2000
      ncnods=nunods
c
c  assign program nodes
c
  200 do 280 id=1,50
      if (nnods(id).eq.0) go to 280
      loc=locate(id)
  210 if (loc.eq.0) go to 280
      if (id.le.4) go to 220
      if (id.le.8) go to 250
      if (id.eq.19) go to 265
      if (id.le.40) go to 220
      if (id.le.43) go to 240
  220 jstop=loc+nnods(id)-1
      do 230 j=loc,jstop
      call getnod(nodplc(j+2))
  230 continue
      go to 270
  240 if (nodplc(loc+5).eq.0) go to 220
      go to 270
  250 call getnod(nodplc(loc+2))
      call getnod(nodplc(loc+3))
      if (id.ge.7) go to 270
      locp=nodplc(loc+id+1)
      nssnod=2*nodplc(loc+4)
  255 do 260 j=1,nssnod
      call getnod(nodplc(locp+j))
  260 continue
      go to 270
  265 locp=nodplc(loc+2)
      call sizmem(nodplc(loc+2),nssnod)
      go to 255
  270 loc=nodplc(loc)
      go to 210
  280 continue
c
c  check and set .nodeset nodes to their internal values
c
      call sizmem(nsnod,nic)
      if(nic.eq.0) go to 300
      do 290 i=1,nic
      call getnod(nodplc(nsnod+i))
  290 continue
c
c   check and set .ic nodes to their internal values
c
  300 call sizmem(icnod,nic)
      if(nic.eq.0) go to 320
      do 310 i=1,nic
      call getnod(nodplc(icnod+i))
  310 continue
  320 if (nogo.ne.0) go to 2000
c
c  expand subcircuit calls
c
      call subckt
      if (nogo.ne.0) go to 2000
      if (ncnods.ge.2) go to 400
      write (iofile,321)
  321 format('0*error*:  circuit has no nodes'/)
      nogo=1
      go to 2000
  400 numnod=ncnods
c
c  link unsatisfied references
c
      call lnkref
      if (nogo.ne.0) go to 2000
c
c  generate subcircuit element names
c
      if (jelcnt(19).eq.0) go to 530
      do 520 id=1,24
      loc=locate(id)
  510 if (loc.eq.0) go to 520
      call subnam(loc)
      loc=nodplc(loc)
      go to 510
  520 continue
c
c  translate node initial conditions to device initial conditions
c  (capacitance, diode, bjt, jfet and mosfet only) when uic is
c  specified on the .tran card
c
  530 if (nosolv.le.0) go to 600
      call sizmem(icnod,nic)
      if(nic.eq.0) go to 600
      call getm8(lvnim1,numnod)
      call zero8(value(lvnim1+1),numnod)
      do 535 i=1,nic
      node=nodplc(icnod+i)
  535 value(lvnim1+node)=value(icval+i)
      loc=locate(2)
  540 if(loc.eq.0) go to 550
      locv=nodplc(loc+1)
      if(value(locv+2).ne.0.0d0) go to 545
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      value(locv+2)=value(lvnim1+node1)-value(lvnim1+node2)
  545 loc=nodplc(loc)
      go to 540
  550 loc=locate(11)
  555 if(loc.eq.0) go to 565
      locv=nodplc(loc+1)
      if(value(locv+2).ne.0.0d0) go to 560
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      value(locv+2)=value(lvnim1+node1)-value(lvnim1+node2)
  560 loc=nodplc(loc)
      go to 555
  565 loc=locate(12)
  570 if(loc.eq.0) go to 580
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      if(value(locv+2).eq.0.0d0) value(locv+2)=value(lvnim1+node2)-
     1  value(lvnim1+node3)
      if(value(locv+3).eq.0.0d0) value(locv+3)=value(lvnim1+node1)-
     1  value(lvnim1+node3)
      loc=nodplc(loc)
      go to 570
  580 loc=locate(13)
  585 if(loc.eq.0) go to 590
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      if(value(locv+2).eq.0.0d0) value(locv+2)=value(lvnim1+node1)-
     1  value(lvnim1+node3)
      if(value(locv+3).eq.0.0d0) value(locv+3)=value(lvnim1+node2)-
     1  value(lvnim1+node3)
      loc=nodplc(loc)
      go to 585
  590 loc=locate(14)
  595 if(loc.eq.0) go to 598
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      if(value(locv+5).eq.0.0d0) value(locv+5)=value(lvnim1+node1)-
     1  value(lvnim1+node3)
      if(value(locv+6).eq.0.0d0) value(locv+6)=value(lvnim1+node2)-
     1  value(lvnim1+node3)
      if(value(locv+7).eq.0.0d0) value(locv+7)=value(lvnim1+node4)-
     1  value(lvnim1+node3)
      loc=nodplc(loc)
      go to 595
  598 call clrmem(lvnim1)
c
c  process sources
c
  600 if (jtrflg.eq.0) go to 700
      do 690 id=9,10
      loc=locate(id)
  610 if (loc.eq.0) go to 690
      locv=nodplc(loc+1)
      locp=nodplc(loc+5)
      jtype=nodplc(loc+4)+1
      go to (680,620,630,640,650,675), jtype
  620 value(locp+3)=dmax1(value(locp+3),0.0d0)
      if (value(locp+4).le.0.0d0) value(locp+4)=tstep
      if (value(locp+5).le.0.0d0) value(locp+5)=tstep
      if (value(locp+6).le.0.0d0) value(locp+6)=tstop
      if (value(locp+7).le.0.0d0) value(locp+7)=tstop
      temp=value(locp+4)+value(locp+5)+value(locp+6)
      value(locp+7)=dmax1(value(locp+7),temp)
      value(locv+1)=value(locp+1)
      go to 680
  630 if (value(locp+3).le.0.0d0) value(locp+3)=1.0d0/tstop
      value(locp+4)=dmax1(value(locp+4),0.0d0)
      value(locv+1)=value(locp+1)
      go to 680
  640 value(locp+3)=dmax1(value(locp+3),0.0d0)
      if (value(locp+4).le.0.0d0) value(locp+4)=tstep
      if (value(locp+5).le.value(locp+3))
     1   value(locp+5)=value(locp+3)+tstep
      if (value(locp+6).le.0.0d0) value(locp+6)=tstep
      value(locv+1)=value(locp+1)
      go to 680
  650 value(locp+1)=dmin1(dmax1(value(locp+1),0.0d0),tstop)
      iknt=1
      call sizmem(nodplc(loc+5),nump)
  660 temp=value(locp+iknt)
      if (value(locp+iknt+2).eq.0.0d0) go to 670
      if (value(locp+iknt+2).ge.tstop) go to 670
      value(locp+iknt+2)=dmax1(value(locp+iknt+2),temp)
      if(temp.ne.value(locp+iknt+2)) go to 665
      write(iofile,661) value(locv)
  661 format('0*error*:  element ',a8,' piecewise linear source table no
     1t increasing in time')
      nogo=1
  665 iknt=iknt+2
      if (iknt.lt.nump) go to 660
  670 value(locp+iknt+2)=tstop
      value(locv+1)=value(locp+2)
      call relmem(nodplc(loc+5),nump-iknt-3)
      go to 680
  675 if (value(locp+3).le.0.0d0) value(locp+3)=1.0d0/tstop
      if (value(locp+5).le.0.0d0) value(locp+5)=1.0d0/tstop
      value(locv+1)=value(locp+1)
  680 loc=nodplc(loc)
      go to 610
  690 continue
c
c  use default values for mos device geometries if not specified
c
  700 loc=locate(14)
  710 if(loc.eq.0) go to 720
      locv=nodplc(loc+1)
      if(value(locv+1).le.0.0d0) value(locv+1)=defl
      if(value(locv+2).le.0.0d0) value(locv+2)=defw
      if(value(locv+3).le.0.0d0) value(locv+3)=defad
      if(value(locv+4).le.0.0d0) value(locv+4)=defas
      loc=nodplc(loc)
      go to 710
c
c  print listing of elements, process device models,
c  and check topology
c
  720 if (iprntl.eq.0) go to 730
      call elprnt
  730 call topchk
      call modchk
      if (nogo.ne.0) go to 2000
c
c  invert resistance values
c
  800 loc=locate(1)
  810 if (loc.eq.0) go to 900
      locv=nodplc(loc+1)
      value(locv+1)=1.0d0/value(locv+2)
      loc=nodplc(loc)
      go to 810
c
c  process mutual inductors
c
  900 loc=locate(4)
  910 if (loc.eq.0) go to 1000
      locv=nodplc(loc+1)
      nl1=nodplc(loc+2)
      lptr1=nodplc(nl1+1)
      nl2=nodplc(loc+3)
      lptr2=nodplc(nl2+1)
      value(locv+1)=value(locv+1)*dsqrt(value(lptr1+1)*value(lptr2+1))
      loc=nodplc(loc)
      go to 910
c
c  limit delmax  if transmission lines in circuit
c
 1000 if (jtrflg.eq.0) go to 1200
      tdmax=0.0d0
      loc=locate(17)
 1010 if (loc.eq.0) go to 1200
      locv=nodplc(loc+1)
      delmax=dmin1(delmax,value(locv+2)/2.0d0)
      tdmax=dmax1(tdmax,value(locv+2))
      loc=nodplc(loc)
      go to 1010
c
c  process source parameters
c
 1200 numbkp=0
      if (jtrflg.eq.0) go to 1205
      tol=1.0d-2*delmax
      numbkp=2
      call getm8(lsbkpt,numbkp)
      value(lsbkpt+1)=0.0d0
      value(lsbkpt+2)=tstop
 1205 do 1290 id=9,10
      loc=locate(id)
 1210 if (loc.eq.0) go to 1290
      locv=nodplc(loc+1)
      locp=nodplc(loc+5)
      temp=value(locv+3)/rad
      value(locv+3)=value(locv+2)*dsin(temp)
      value(locv+2)=value(locv+2)*dcos(temp)
      if (jtrflg.eq.0) go to 1280
      jtype=nodplc(loc+4)+1
      go to (1280,1220,1230,1235,1240,1260), jtype
 1220 value(locp+4)=value(locp+4)+value(locp+3)
      temp=value(locp+5)
      value(locp+5)=value(locp+4)+value(locp+6)
      value(locp+6)=value(locp+5)+temp
      time=0.0d0
 1225 call extmem(lsbkpt,4)
      value(lsbkpt+numbkp+1)=value(locp+3)+time
      value(lsbkpt+numbkp+2)=value(locp+4)+time
      value(lsbkpt+numbkp+3)=value(locp+5)+time
      value(lsbkpt+numbkp+4)=value(locp+6)+time
      numbkp=numbkp+4
      time=time+value(locp+7)
      if (time.ge.tstop) go to 1280
      go to 1225
 1230 value(locp+3)=value(locp+3)*twopi
      call extmem(lsbkpt,1)
 1231 value(lsbkpt+numbkp+1)=value(locp+4)
      numbkp=numbkp+1
      go to 1280
 1235 call extmem(lsbkpt,2)
      value(lsbkpt+numbkp+1)=value(locp+3)
      value(lsbkpt+numbkp+2)=value(locp+5)
      numbkp=numbkp+2
      go to 1280
 1240 iknt=1
      call sizmem(nodplc(loc+5),nump)
 1250 call extmem(lsbkpt,1)
      value(lsbkpt+numbkp+1)=value(locp+iknt)
      numbkp=numbkp+1
      iknt=iknt+2
      if (iknt.le.nump) go to 1250
      go to 1280
 1260 value(locp+3)=value(locp+3)*twopi
      value(locp+5)=value(locp+5)*twopi
 1280 loc=nodplc(loc)
      go to 1210
 1290 continue
c
c  augment breakpoint table for transmission line delays
c
      if (jtrflg.eq.0) go to 1300
      loc=locate(17)
 1292 if (loc.eq.0) go to 1300
      locv=nodplc(loc+1)
      td=value(locv+2)
      ntemp=numbkp
      do 1296 ibkp=1,ntemp
      time=value(lsbkpt+ibkp)
 1294 time=time+td
      if (time.ge.tstop) go to 1296
      call extmem(lsbkpt,1)
      value(lsbkpt+numbkp+1)=time
      numbkp=numbkp+1
      go to 1294
 1296 continue
      call shlsrt(value(lsbkpt+1),numbkp)
      nbkpt=1
      do 1298 i=2,numbkp
      if ((value(lsbkpt+i)-value(lsbkpt+nbkpt)).lt.tol) go to 1298
      nbkpt=nbkpt+1
      value(lsbkpt+nbkpt)=value(lsbkpt+i)
      if (value(lsbkpt+nbkpt).ge.tstop) go to 1299
 1298 continue
 1299 call relmem(lsbkpt,numbkp-nbkpt)
      numbkp=nbkpt
      value(lsbkpt+numbkp)=dmax1(value(lsbkpt+numbkp),tstop)
      loc=nodplc(loc)
      go to 1292
c
c  finish breakpoint table
c
 1300 if (jtrflg.eq.0) go to 1600
      call extmem(lsbkpt,1)
      value(lsbkpt+numbkp+1)=tstop
      numbkp=numbkp+1
      call shlsrt(value(lsbkpt+1),numbkp)
      nbkpt=1
      do 1310 i=2,numbkp
      if ((value(lsbkpt+i)-value(lsbkpt+nbkpt)).lt.tol) go to 1310
      nbkpt=nbkpt+1
      value(lsbkpt+nbkpt)=value(lsbkpt+i)
      if (value(lsbkpt+nbkpt).ge.tstop) go to 1320
 1310 continue
 1320 call relmem(lsbkpt,numbkp-nbkpt)
      numbkp=nbkpt
      value(lsbkpt+numbkp)=dmax1(value(lsbkpt+numbkp),tstop)
c
c  print option summary
c
 1600 if (iprnto.eq.0) go to 1700
      call title(0,lwidth,1,titlop)
      write (iofile,1601) gmin,reltol,abstol,vntol,lvlcod,itl1,itl2
 1601 format('0dc analysis -',/,
     1   '0    gmin   = ',1pd10.3,/,
     2   '     reltol = ',  d10.3,/,
     3   '     abstol = ',  d10.3,/,
     4   '     vntol  = ',  d10.3,/,
     5   '     lvlcod = ',     i6,/,
     6   '     itl1   = ',     i6,/,
     7   '     itl2   = ',     i6,/)
      write (iofile,1605) pivtol,pivrel
 1605 format(
     1   '     pivtol = ',1pd10.3,/,
     2   '     pivrel = ',  d10.3)
      write (iofile,1611) aname(method),maxord,chgtol,trtol,lvltim,xmu,
     1   itl3,itl4,itl5
 1611 format('0transient analysis -',/,
     1   '0    method =  ',a8,/,
     2   '     maxord = ',     i6,/,
     3   '     chgtol = ',1pd10.3,/,
     4   '     trtol  = ',  d10.3,/,
     5   '     lvltim = ',     i6,/,
     6   '     mu     = ',0pf10.3,/,
     7   '     itl3   = ',     i6,/,
     8   '     itl4   = ',     i6,/,
     9   '     itl5   = ',     i6,/)
      write (iofile,1621) limpts,limtim,maxtim,numdgt,value(itemps+1),
     1   defl,defw,defad,defas
 1621 format('0miscellaneous -',/,
     1   '0    limpts = ',     i6,/,
     2   '     limtim = ',     i6,/,
     3   '     cptime = ',     i9,/,
     4   '     numdgt = ',     i6,/,
     5   '     tnom   = ',0pf10.3,/,
     6   '     defl   = ',1pd10.3,/,
     7   '     defw   = ',d10.3,/,
     8   '     defad  = ',d10.3,/,
     9   '     defas  = ',d10.3)
c
c  miscellaneous error checking
c
 1700 if (icvflg.eq.0) go to 1720
      if (icvflg.le.limpts) go to 1710
      icvflg=0
      write (iofile,1701) limpts,alsdc
 1701 format('0warning:  more than ',i5,' points for ',a4,' analysis,',/
     11x,'analysis omitted.  this limit may be overridden using the ',/
     21x,'limpts parameter on the .option card'/)
      go to 1720
 1710 if ((jelcnt(31)+jelcnt(36)).gt.0) go to 1720
      if(ipostp.ne.0) go to 1720
      icvflg=0
      write (iofile,1711) alsdc
 1711 format('0warning:  no ',a4,' outputs specified .',
     1  '.. analysis omitted'/)
 1720 if (jtrflg.eq.0) go to 1740
      if (method.eq.1) maxord=2
      if ((method.eq.2).and.(maxord.ge.3)) lvltim=2
      if (jtrflg.le.limpts) go to 1730
      jtrflg=0
      write (iofile,1701) limpts,alstr
      go to 1740
 1730 if ((jelcnt(32)+jelcnt(37)+nfour).gt.0) go to 1735
      if(ipostp.ne.0) go to 1735
      jtrflg=0
      write (iofile,1711) alstr
      go to 1740
 1735 if (nfour.eq.0) go to 1740
      forprd=1.0d0/forfre
      if ((tstop-forprd).ge.(tstart-1.0d-12)) go to 1740
      nfour=0
      call clrmem(ifour)
      write (iofile,1736)
 1736 format('0warning:  fourier analysis fundamental frequency is incom
     1patible with'/11x,'transient analysis print interval ... fourier a
     2nalysis omitted'/)
 1740 if (jacflg.eq.0) go to 1800
      if (jacflg.le.limpts) go to 1750
      jacflg=0
      write (iofile,1701) limpts,alsac
      go to 1800
 1750 if ((jelcnt(33)+jelcnt(34)+jelcnt(35)+jelcnt(38)+jelcnt(39)
     1   +jelcnt(40)+idist+inoise).gt.0) go to 1800
      if(ipostp.ne.0) go to 1800
      jacflg=0
      write (iofile,1711) alsac
c
c  sequence through the output lists
c
 1800 do 1820 id=41,45
      if (id.le.43) numout=1
      loc=locate(id)
 1810 if (loc.eq.0) go to 1820
      numout=numout+1
      nodplc(loc+4)=numout
      loc=nodplc(loc)
      go to 1810
 1820 continue
c
c   increase number of .prints if too many outputs for output line-width
c
      ifwdth=max0(numdgt-1,0)+9
      noprln=min0(8,(lwidth-12)/ifwdth)
      do 1860 id=31,35
      loc=locate(id)
 1830 if(loc.eq.0) go to 1860
      noprex=nodplc(loc+3)-noprln
      if(noprex.le.0) go to 1850
      nodplc(loc+3)=noprln
      call find(dble(jelcnt(id)),id,locnew,1)
      nodplc(locnew+2)=nodplc(loc+2)
      nodplc(locnew+3)=noprex
      call copy4(nodplc(loc+2*noprln+4),nodplc(locnew+4),2*noprex)
 1850 loc=nodplc(loc)
      go to 1830
 1860 continue
c
c  exit
c
 2000 call second(t2)
      rstats(1)=rstats(1)+t2-t1
      return
      end
