      subroutine card
      implicit double precision (a-h,o-z)
c
c     this routine scans the input lines, storing each field into the
c tables ifield, idelim, icolum, and icode.  with the exception of the
c '.end' line, card always reads the next line to check for a possible
c continuation before it exits.
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
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      dimension adigit(10)
      data adigit / 1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9 /
      data ablnk,aper,aplus,aminus,astk / 1h , 1h., 1h+, 1h-, 1h* /
      data bg,ak,au,an,ap,ae,am,af,at /1hg,1hk,1hu,1hn,1hp,1he,1hm,
     1  1hf,1ht/
      data ai / 1hi /
      data alprn, arprn, aequal / 1h(, 1h), 1h= /
      data aend / 4h.end /
c
c      note:  the value of the function *nxtchr* (used extensively in
c this routine) is as follows:
c
c                    <0:  end-of-line
c                    =0:  delimiter found
c                    >0:  non-delimiter found
c
      numfld=0
      nofld=10
      go to 20
c
c  read next card
c
   10 nofld=10
      call getlin
      if (keof.eq.0) go to 20
c...  error:  unexpected end-of-file condition on input
   15 keof=1
      nofld=1
      numfld=0
      igoof=1
      write (iofile,16)
   16 format('0*error*:  .end card missing'/)
      go to 1000
c
c  eliminate trailing blanks rapidly
c
   20 if (afield(nofld).ne.ablnk) go to 40
      if (nofld.eq.1) go to 30
      nofld=nofld-1
      go to 20
c...  write blank card
   30 write (iofile,31)
   31 format(1x)
      go to 10
c...  copy the card to output listing
   40 write (iofile,41) (afield(i),i=1,nofld)
   41 format(1x,10a8)
c
c  initialization for new card
c
   45 kntrc=0
      kntlim=min0(8*nofld,iwidth)
c
c  fetch first non-delimiter (see routine *nxtchr* for list)
c
   50 if (nxtchr(0)) 600,50,60
c...  check for comment (leading asterisk)
   60 if (achar.eq.astk) go to 10
      go to 100
c
c  fetch next character
c
   70 if (nxtchr(0)) 600,80,100
c
c  two consecutive delimiters imply numeric zero unless the delimiter
c  is a blank or parenthesis.
c
   80 if (achar.eq.ablnk) go to 70
      if (achar.eq.alprn) go to 70
      if (achar.eq.arprn) go to 70
      if (achar.eq.aequal) go to 70
c...  check for sufficient space in storage arrays
      if (numfld.lt.insize-1) go to 90
      call extmem(ifield,50)
      call extmem(icode,50)
      call extmem(idelim,50)
      call extmem(icolum,50)
      insize=insize+50
   90 numfld=numfld+1
      value(ifield+numfld)=0.0d0
      nodplc(icode+numfld)=0
      value(idelim+numfld)=achar
      nodplc(icolum+numfld)=kntrc
      go to 70
c
c  check for sufficient space in storage arrays
c
  100 if (numfld.lt.insize-1) go to 110
      call extmem(ifield,50)
      call extmem(icode,50)
      call extmem(idelim,50)
      call extmem(icolum,50)
      insize=insize+50
c
c  begin scan of next field
c
c...  initialization
  110 jdelim=0
      xsign=1.0d0
      xmant=0.0d0
      idec=0
      iexp=0
c...  check for leading plus or minus sign
      if (achar.eq.aplus) go to 210
      if (achar.eq.aminus) go to 200
c...  finish initialization
      anam=ablnk
      kchr=1
c...  an isolated period indicates that a continuation card follows
      if (achar.ne.aper) go to 120
c...  alter initialization slightly if leading period found
      idec=1
      iexp=-1
      anam=aper
      kchr=2
c...  now take a look at the next character
      if (nxtchr(0)) 10,10,120
c
c  test for number (any digit)
c
  120 do 130 i=1,10
      if (achar.ne.adigit(i)) go to 130
      xmant=dble(i-1)
      go to 210
  130 continue
c
c  assemble name
c
      numfld=numfld+1
      call move(anam,kchr,achar,1,1)
      kchr=kchr+1
      do 150 i=kchr,8
      if (nxtchr(0)) 160,160,140
  140 call move(anam,i,achar,1,1)
  150 continue
      go to 170
  160 jdelim=1
  170 value(ifield+numfld)=anam
      nodplc(icode+numfld)=1
      nodplc(icolum+numfld)=kntrc
c...  no '+' format continuation possible for .end card
      if (numfld.ge.2) go to 400
      if (anam.ne.aend) go to 400
      nodplc(icode+numfld+1)=-1
      go to 1000
c
c  process number
c
c...  take note of leading minus sign
  200 xsign=-1.0d0
c...  take a look at the next character
  210 if (nxtchr(0)) 335,335,220
c...  test for digit
  220 do 230 i=1,10
      if (achar.ne.adigit(i)) go to 230
      xmant=xmant*10.0d0+dble(i-1)
      if (idec.eq.0) go to 210
      iexp=iexp-1
      go to 210
  230 continue
c
c  check for decimal point
c
      if (achar.ne.aper) go to 240
c...  make certain that this is the first one found
      if (idec.ne.0) go to 500
      idec=1
      go to 210
c
c  test for exponent
c
  240 if (achar.ne.ae) go to 300
      if (nxtchr(0)) 335,335,250
  250 itemp=0
      isign=1
c...  check for possible leading sign on exponent
      if (achar.eq.aplus) go to 260
      if (achar.ne.aminus) go to 270
      isign=-1
  260 if (nxtchr(0)) 285,285,270
c...  test for digit
  270 do 280 i=1,10
      if (achar.ne.adigit(i)) go to 280
      itemp=itemp*10+i-1
      go to 260
  280 continue
      go to 290
  285 jdelim=1
c...  correct internal exponent
  290 iexp=iexp+isign*itemp
      go to 340
c
c  test for scale factor
c
  300 if (achar.ne.am) go to 330
c...  special check for *me* (as distinguished from *m*)
      if (nxtchr(0)) 320,320,310
  310 if (achar.ne.ae) go to 315
      iexp=iexp+6
      go to 340
  315 if (achar.ne.ai) go to 325
      xmant=xmant*25.4d-6
      go to 340
  320 jdelim=1
  325 iexp=iexp-3
      go to 340
  330 if (achar.eq.at) iexp=iexp+12
      if (achar.eq.bg) iexp=iexp+9
      if (achar.eq.ak) iexp=iexp+3
      if (achar.eq.au) iexp=iexp-6
      if (achar.eq.an) iexp=iexp-9
      if (achar.eq.ap) iexp=iexp-12
      if (achar.eq.af) iexp=iexp-15
      go to 340
  335 jdelim=1
c
c  assemble the final number
c
  340 if (xmant.eq.0.0d0) go to 350
      if (iexp.eq.0) go to 350
      if (iabs(iexp).ge.201) go to 500
      xmant=xmant*dexp(dble(iexp)*xlog10)
      if (xmant.gt.1.0d+35) go to 500
      if (xmant.lt.1.0d-35) go to 500
  350 numfld=numfld+1
      value(ifield+numfld)=dsign(xmant,xsign)
      nodplc(icode+numfld)=0
      nodplc(icolum+numfld)=kntrc
c
c  skip to non-blank delimiter (if necessary)
c
  400 if (jdelim.eq.0) go to 440
  410 value(idelim+numfld)=achar
      if (achar.ne.ablnk) go to 70
      if (nxtchr(0)) 450,410,420
  420 kntrc=kntrc-1
      go to 70
  440 if (nxtchr(0)) 450,410,440
  450 value(idelim+numfld)=achar
      go to 600
c
c  errors
c
  500 write (iofile,501) kntrc
  501 format('0*error*:  illegal number -- scan stopped at column ',i3/)
      igoof=1
      numfld=numfld+1
      value(ifield+numfld)=0.0d0
      nodplc(icode+numfld)=0
      value(idelim+numfld)=achar
      nodplc(icolum+numfld)=kntrc
c
c  finished
c
  600 nodplc(icode+numfld+1)=-1
c
c  check next line for possible continuation
c
  610 call getlin
      if (keof.eq.1) go to 15
      nofld=10
  620 if (afield(nofld).ne.ablnk) go to 630
      if (nofld.eq.1) go to 650
      nofld=nofld-1
      go to 620
  630 kntrc=0
      kntlim=min0(8*nofld,iwidth)
c...  continuation line has a '+' as first non-delimiter on card
  632 if(nxtchr(0)) 650,632,634
  634 if(achar.ne.aplus) go to 640
      write(iofile,41) (afield(i),i=1,nofld)
      go to 70
  640 if (achar.ne.astk) go to 1000
  650 write (iofile,41) (afield(i),i=1,nofld)
      go to 610
 1000 return
      end
