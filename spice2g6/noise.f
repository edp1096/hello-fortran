      subroutine noise(loco)
      implicit double precision (a-h,o-z)
c
c     this routine computes the noise due to various circuit elements.
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
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension vno1(12),vno2(12),vno3(12),vno4(12),vno5(12),vno6(12)
      dimension vntot(12),anam(12),string(5)
      real v,vreal,vimag
      dimension titln(4),v(2)
      dimension afmt1(3),afmt2(3)
      complex cval,c(1)
      equivalence (c(1),v(1),cval)
      equivalence (v(1),vreal),(v(2),vimag)
      data titln / 8hnoise an, 8halysis  , 8h        , 8h         /
      data alsrb,alsrc,alsre,alsrs,alsrd / 2hrb,2hrc,2hre,2hrs,2hrd /
      data alsib,alsic,alsid,alsfn / 2hib,2hic,2hid,2hfn /
      data alstot / 5htotal /
      data aslash,ablnk / 1h/, 1h  /
      data afmt1 /8h(////,11,8hx,  (2x,,8ha8))    /
      data afmt2 /8h(1h0,a8,,8h1p  d10.,8h3)      /
c
c
c.. fix-up formats
      kntr=12
      if(lwidth.le.80) kntr=7
      ipos=11
      call move(afmt1,ipos,ablnk,1,2)
      call alfnum(kntr,afmt1,ipos)
      ipos=11
      call move(afmt2,ipos,ablnk,1,2)
      call alfnum(kntr,afmt2,ipos)
      nprnt=0
      freq=omega/twopi
      if (icalc.ge.2) go to 10
      fourkt=4.0d0*charge*vt
      twoq=2.0d0*charge
      noposo=nodplc(nosout+2)
      nonego=nodplc(nosout+3)
      kntlim=lwidth/11
      nkntr=1
   10 if (nosprt.eq.0) go to 30
      if (nkntr.gt.icalc) go to 30
      nprnt=1
      nkntr=nkntr+nosprt
      call title(0,lwidth,1,titln)
      write (iofile,16) freq
   16 format('0    frequency = ',1pd10.3,' hz'/)
c
c  obtain adjoint circuit solution
c
   30 vnrms=0.0d0
      cval=cvalue(lcvn+noposo)-cvalue(lcvn+nonego)
      vout=dsqrt(dble(vreal*vreal)+dble(vimag*vimag))
      vout=dmax1(vout,1.0d-20)
      call zero8(value(lvn+1),nstop)
      call zero8(value(imvn+1),nstop)
      value(lvn+noposo)=-1.0d0
      value(lvn+nonego)=+1.0d0
      call acasol
c
c  resistors
c
      if (jelcnt(1).eq.0) go to 200
      ititle=0
   91 format(//'0**** resistor squared noise voltages (sq v/hz)')
  100 loc=locate(1)
      kntr=0
  110 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 130
      kntr=kntr+1
      locv=nodplc(loc+1)
      anam(kntr)=value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      cval=cvalue(lcvn+node1)-cvalue(lcvn+node2)
      vntot(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locv+1)
      vnrms=vnrms+vntot(kntr)
      if (kntr.ge.kntlim) go to 140
  120 loc=nodplc(loc)
      go to 110
  130 if (kntr.eq.0) go to 200
  140 if (nprnt.eq.0) go to 160
      if (ititle.eq.0) write (iofile,91)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) alstot,(vntot(i),i=1,kntr)
  160 kntr=0
      if (loc.ne.0) go to 120
c
c  diodes
c
  200 if (jelcnt(11).eq.0) go to 300
      ititle=0
  201 format(//'0**** diode squared noise voltages (sq v/hz)')
  210 loc=locate(11)
      kntr=0
  220 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 240
      kntr=kntr+1
      locv=nodplc(loc+1)
      anam(kntr)=value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      loct=nodplc(loc+11)
      area=value(locv+1)
      fnk=value(locm+10)
      fna=value(locm+11)
c
c  ohmic resistance
c
      cval=cvalue(lcvn+node1)-cvalue(lcvn+node3)
      vno1(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locm+2)*area
c
c  junction shot noise and flicker noise
c
      cval=cvalue(lcvn+node3)-cvalue(lcvn+node2)
      vtemp=dble(vreal*vreal)+dble(vimag*vimag)
      arg=dmax1(dabs(value(lx0+loct+1)),1.0d-20)
      vno2(kntr)=vtemp*twoq*arg
      vno3(kntr)=vtemp*fnk*dexp(fna*dlog(arg))/freq
      vntot(kntr)=vno1(kntr)+vno2(kntr)+vno3(kntr)
      vnrms=vnrms+vntot(kntr)
      if (kntr.ge.kntlim) go to 250
  230 loc=nodplc(loc)
      go to 220
  240 if (kntr.eq.0) go to 300
  250 if (nprnt.eq.0) go to 260
      if (ititle.eq.0) write (iofile,201)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) alsrs,(vno1(i),i=1,kntr)
      write (iofile,afmt2) alsid,(vno2(i),i=1,kntr)
      write (iofile,afmt2) alsfn,(vno3(i),i=1,kntr)
      write (iofile,afmt2) alstot,(vntot(i),i=1,kntr)
  260 kntr=0
      if (loc.ne.0) go to 230
c
c  bipolar junction transistors
c
  300 if (jelcnt(12).eq.0) go to 400
      ititle=0
  301 format(//'0**** transistor squared noise voltages (sq v/hz)')
  310 loc=locate(12)
      kntr=0
  320 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 340
      kntr=kntr+1
      locv=nodplc(loc+1)
      anam(kntr)=value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      loct=nodplc(loc+22)
      area=value(locv+1)
      fnk=value(locm+44)
      fna=value(locm+45)
c
c  extrinsic resistances
c
c...  base resistance
      cval=cvalue(lcvn+node2)-cvalue(lcvn+node5)
      vno1(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1  *fourkt*value(lx0+loct+16)
c...  collector resistance
      cval=cvalue(lcvn+node1)-cvalue(lcvn+node4)
      vno2(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locm+20)*area
c...  emitter resistance
      cval=cvalue(lcvn+node3)-cvalue(lcvn+node6)
      vno3(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locm+19)*area
c
c  base current shot noise and flicker noise
c
      cval=cvalue(lcvn+node5)-cvalue(lcvn+node6)
      vtemp=dble(vreal*vreal)+dble(vimag*vimag)
      arg=dmax1(dabs(value(lx0+loct+3)),1.0d-20)
      vno4(kntr)=vtemp*twoq*arg
      vno5(kntr)=vtemp*fnk*dexp(fna*dlog(arg))/freq
c
c  collector current shot noise
c
      cval=cvalue(lcvn+node4)-cvalue(lcvn+node6)
      vno6(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *twoq*dabs(value(lx0+loct+2))
      vntot(kntr)=vno1(kntr)+vno2(kntr)+vno3(kntr)+vno4(kntr)+vno5(kntr)
     1   +vno6(kntr)
      vnrms=vnrms+vntot(kntr)
      if (kntr.ge.kntlim) go to 350
  330 loc=nodplc(loc)
      go to 320
  340 if (kntr.eq.0) go to 400
  350 if (nprnt.eq.0) go to 360
      if (ititle.eq.0) write (iofile,301)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) alsrb,(vno1(i),i=1,kntr)
      write (iofile,afmt2) alsrc,(vno2(i),i=1,kntr)
      write (iofile,afmt2) alsre,(vno3(i),i=1,kntr)
      write (iofile,afmt2) alsib,(vno4(i),i=1,kntr)
      write (iofile,afmt2) alsic,(vno6(i),i=1,kntr)
      write (iofile,afmt2) alsfn,(vno5(i),i=1,kntr)
      write (iofile,afmt2) alstot,(vntot(i),i=1,kntr)
  360 kntr=0
      if (loc.ne.0) go to 330
c
c  jfets
c
  400 if (jelcnt(13).eq.0) go to 500
      ititle=0
  401 format(//'0**** jfet squared noise voltages (sq v/hz)')
  410 loc=locate(13)
      kntr=0
  420 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 440
      kntr=kntr+1
      locv=nodplc(loc+1)
      anam(kntr)=value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      locm=nodplc(loc+7)
      locm=nodplc(locm+1)
      loct=nodplc(loc+19)
      area=value(locv+1)
      fnk=value(locm+10)
      fna=value(locm+11)
c
c  extrinsic resistances
c
c...  drain resistance
      cval=cvalue(lcvn+node1)-cvalue(lcvn+node4)
      vno1(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locm+4)*area
c...  source resistance
      cval=cvalue(lcvn+node3)-cvalue(lcvn+node5)
      vno2(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))
     1   *fourkt*value(locm+5)*area
c
c  drain current shot noise and flicker noise
c
      cval=cvalue(lcvn+node4)-cvalue(lcvn+node5)
      vtemp=dble(vreal*vreal)+dble(vimag*vimag)
      vno3(kntr)=vtemp*fourkt*2.0d0*dabs(value(lx0+loct+5))/3.0d0
      arg=dmax1(dabs(value(lx0+loct+3)),1.0d-20)
      vno4(kntr)=vtemp*fnk*dexp(fna*dlog(arg))/freq
      vntot(kntr)=vno1(kntr)+vno2(kntr)+vno3(kntr)+vno4(kntr)
      vnrms=vnrms+vntot(kntr)
      if (kntr.ge.kntlim) go to 450
  430 loc=nodplc(loc)
      go to 420
  440 if (kntr.eq.0) go to 500
  450 if (nprnt.eq.0) go to 460
      if (ititle.eq.0) write (iofile,401)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) alsrd,(vno1(i),i=1,kntr)
      write (iofile,afmt2) alsrs,(vno2(i),i=1,kntr)
      write (iofile,afmt2) alsid,(vno3(i),i=1,kntr)
      write (iofile,afmt2) alsfn,(vno4(i),i=1,kntr)
      write (iofile,afmt2) alstot,(vntot(i),i=1,kntr)
  460 kntr=0
      if (loc.ne.0) go to 430
c
c  mosfets
c
  500 if (jelcnt(14).eq.0) go to 600
      ititle=0
  501 format(//'0**** mosfet squared noise voltages (sq v/hz)')
  510 loc=locate(14)
      kntr=0
  520 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 540
      kntr=kntr+1
      locv=nodplc(loc+1)
      anam(kntr)=value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      locm=nodplc(loc+8)
      itype=nodplc(locm+2)
      loct=nodplc(loc+26)
      locm=nodplc(locm+1)
      xl=value(locv+1)-2.0d0*value(locm+28)
      xw=value(locv+2)
      cox=value(locm+22)
      if (cox.le.0.0d0) cox=epsox/1.0d-7
      fnk=value(locm+36)
      fna=value(locm+37)
c
c  extrinsic resistances
c
      if ((value(locm+7).le.0.0d0).and.
     1   (value(locm+8).le.0.0d0)) go to 522
      gdpr=value(locm+7)
      gspr=value(locm+8)
      go to 524
  522 gdpr=value(locm+16)/value(locv+13)
      gspr=value(locm+16)/value(locv+14)
c...  drain resistance
  524 cval=cvalue(lcvn+node1)-cvalue(lcvn+node5)
      vno1(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))*fourkt*gdpr
c...  source resistance
      cval=cvalue(lcvn+node3)-cvalue(lcvn+node6)
      vno2(kntr)=(dble(vreal*vreal)+dble(vimag*vimag))*fourkt*gspr
c
c  drain current shot noise and flicker noise
c
      cval=cvalue(lcvn+node5)-cvalue(lcvn+node6)
      vtemp=dble(vreal*vreal)+dble(vimag*vimag)
      gm=value(lx0+loct+7)
      arg=dmax1(dabs(value(lx0+loct+4)),1.0d-20)
      vno3(kntr)=vtemp*fourkt*dabs(gm)/1.5d0
      vno4(kntr)=vtemp*fnk*dexp(fna*dlog(arg))/(freq*cox*xl*xl)
  525 vntot(kntr)=vno1(kntr)+vno2(kntr)+vno3(kntr)+vno4(kntr)
      vnrms=vnrms+vntot(kntr)
      if (kntr.ge.kntlim) go to 550
  530 loc=nodplc(loc)
      go to 520
  540 if (kntr.eq.0) go to 600
  550 if (nprnt.eq.0) go to 560
      if (ititle.eq.0) write (iofile,501)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) alsrd,(vno1(i),i=1,kntr)
      write (iofile,afmt2) alsrs,(vno2(i),i=1,kntr)
      write (iofile,afmt2) alsid,(vno3(i),i=1,kntr)
      write (iofile,afmt2) alsfn,(vno4(i),i=1,kntr)
      write (iofile,afmt2) alstot,(vntot(i),i=1,kntr)
  560 kntr=0
      if (loc.ne.0) go to 530
c
c  compute equivalent input noise voltage
c
  600 vnout=dsqrt(vnrms)
      vnin=vnout/vout
      if (nprnt.eq.0) go to 620
      do 610 i=1,5
      string(i)=ablnk
  610 continue
      ioutyp=1
      ipos=1
      call outnam(nosout,ioutyp,string,ipos)
      call move(string,ipos,aslash,1,1)
      ipos=ipos+1
      locv=nodplc(nosin+1)
      anam1=value(locv)
      call move(string,ipos,anam1,1,8)
      write (iofile,611) vnrms,vnout,string,vout,anam1,vnin
  611 format(////,
     1   '0**** total output noise voltage',9x,'= ',1pd10.3,' sq v/hz'/,
     2   1h0,40x,'= ',d10.3,' v/rt hz'/,
     3   '0     transfer function value:',/,
     4   1h0,7x,4a8,a1,'= ',d10.3,/,
     5   '0     equivalent input noise at ',a8,' = ',d10.3,' /rt hz')
c
c  save noise outputs
c
  620 loc=locate(44)
  630 if (loc.eq.0) go to 1000
      iseq=nodplc(loc+4)
      if (nodplc(loc+5).ne.2) go to 640
      cvalue(loco+iseq)=vnout
      go to 650
  640 cvalue(loco+iseq)=vnin
  650 loc=nodplc(loc)
      go to 630
c
c  finished
c
 1000 return
      end
