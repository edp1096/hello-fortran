c spice version 2g.6  sccsid=dcop.ma 3/15/83
      subroutine dcop
      implicit double precision (a-h,o-z)
c
c     this routine prints out the operating points of the nonlinear
c circuit elements.
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      logical memptr
c
c
      dimension optitl(4)
      dimension anam(12),av1(12),ai1(12),req(12)
      dimension amod(12),vd(12),cap(12)
      dimension cb(12),cc(12),vbe(12),vbc(12),vce(12),rpi(12),
     1   ro(12),cpi(12),cmu(12),betadc(12),betaac(12),ft(12),
     2   ccs(12),cbx(12),rx(12)
      dimension cg(12),vgs(12),vds(12),gds(12),vbs(12),cbd(12),cbs(12),
     2  cgsov(12),cgdov(12),cgbov(12),vth(12),vdsat(12),cd(12),gm(12),
     3  cggb(12),cgdb(12),cgsb(12),cbgb(12),cbdb(12),cbsb(12),
     4  gmb(12)
      dimension cgs(12),cgd(12),cgb(12),cds(12)
      equivalence(cb(1),cg(1)),(cc(1),vgs(1)),(vbe(1),vds(1)),
     1(vbc(1),gds(1)),(vce(1),vbs(1)),(rpi(1),cbd(1)),
     2(ro(1),cbs(1)),(cpi(1),cgsov(1)),(cmu(1),cgdov(1)),
     3(betadc(1),cgbov(1)),(betaac(1),vth(1)),(ft(1),vdsat(1)),
     4(ccs(1),cd(1)),(cbx(1),cggb(1)),(rx(1),cgdb(1))
      equivalence(vd(1),cg(1)),(cap(1),vgs(1)),(av1(1),vds(1)),
     1  (ai1(1),gds(1)),(req(1),vbs(1))
      equivalence (cgs(1),cggb(1)),(cgd(1),cgdb(1)),(cgb(1),cgsb(1)),
     1  (cds(1),cbgb(1))
      dimension afmt1(3),afmt2(2),afmt3(3),afmt4(3)
      data optitl / 8hoperatin, 8hg point , 8hinformat, 8hion      /
      data av,avd,avbe,avbc,avce,avgs,avds,avbs / 1hv,2hvd,3hvbe,3hvbc,
     1   3hvce,3hvgs,3hvds,3hvbs /
      data acntrv,acntri,asrcv,asrci,atrang,atranr,avgain,aigain /
     1   8hv-contrl, 8hi-contrl, 8hv-source, 8hi-source,
     2   8htrans-g , 8htrans-r , 8hv gain  , 8hi gain   /
      data ai,aid,aib,aic,aig / 1hi,2hid,2hib,2hic,2hig /
      data areq,arpi,aro / 3hreq,3hrpi,2hro /
      data acap,acpi,acmu,acgs,acgd,acbd,acbs / 3hcap,3hcpi,3hcmu,3hcgs,
     1   3hcgd,3hcbd,3hcbs /
      data acgsov,acgdov,acgbov /6hcgsovl,6hcgdovl,6hcgbovl/
      data acggb,acgdb,acgsb,acbgb,acbdb,acbsb /7hdqgdvgb,7hdqgdvdb,
     1  7hdqgdvsb,7hdqbdvgb,7hdqbdvdb,7hdqbdvsb/
      data acgb,acds / 3hcgb,3hcds /
      data avth, avdsat / 3hvth, 5hvdsat /
      data agm,agds / 2hgm,3hgds /
      data agmb / 4hgmb /
      data accs,acbx,arx /3hccs,3hcbx,2hrx/
      data abetad,abetaa / 6hbetadc,6hbetaac /
      data aft / 2hft /
c
      data ablnk /1h /
      data afmt1 /8h(//1h0,1,8h0x,  (2x,8h,a8))   /
      data afmt2 /8h(1h ,a8,,8h  f10.3)/
      data afmt3 /8h(1h ,a8,,8h1p  e10.,8h2)      /
      data afmt4 /8h('0model,8h   ',  (,8h2x,a8)) /
c
c.. fix-up the format statements
c
      kntr=12
      if(lwidth.le.80) kntr=7
      ipos=12
      call move(afmt1,ipos,ablnk,1,2)
      call alfnum(kntr,afmt1,ipos)
      ipos=9
      call move(afmt2,ipos,ablnk,1,2)
      call alfnum(kntr,afmt2,ipos)
      ipos=11
      call move(afmt3,ipos,ablnk,1,2)
      call alfnum(kntr,afmt3,ipos)
      ipos=14
      call move(afmt4,ipos,ablnk,1,2)
      call alfnum(kntr,afmt4,ipos)
c
c  compute voltage source currents and power dissipation
c
      call second(t1)
      if ((mode.eq.1).and.(modedc.eq.2).and.(nosolv.ne.0)) go to 700
      power=0.0d0
      if (jelcnt(9).eq.0) go to 50
      ititle=0
   11 format (////5x,'voltage source currents'//5x,'name',
     1   7x,'current'/)
      loc=locate(9)
   20 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 50
      locv=nodplc(loc+1)
      iptr=nodplc(loc+6)
      creal=value(lvnim1+iptr)
      power=power-creal*value(locv+1)
      if (ititle.eq.0) write (iofile,11)
      ititle=1
      write (iofile,21) value(locv),creal
   21 format (/5x,a8,1x,1pd10.3)
   30 loc=nodplc(loc)
      go to 20
   50 loc=locate(10)
   60 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 90
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      power=power-value(locv+1)
     1   *(value(lvnim1+node1)-value(lvnim1+node2))
      loc=nodplc(loc)
      go to 60
   90 write (iofile,91) power
   91 format (//5x,'total power dissipation  ',1pd9.2,'  watts')
c
c  small signal device parameters
c
      numdev=jelcnt(5)+jelcnt(6)+jelcnt(7)+jelcnt(8)+jelcnt(11)
     1   +jelcnt(12)+jelcnt(13)+jelcnt(14)
      if (numdev.eq.0) go to 600
      call title(0,lwidth,1,optitl)
      kntlim=lwidth/11
c
c  nonlinear voltage controlled current sources
c
      if (jelcnt(5).eq.0) go to 175
      ititle=0
  111 format(1h0,/,'0**** voltage-controlled current sources')
      loc=locate(5)
      kntr=0
  120 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 140
      kntr=kntr+1
      locv=nodplc(loc+1)
      loct=lx0+nodplc(loc+12)
      anam(kntr)=value(locv)
      ai1(kntr)=value(loct)
      if (kntr.ge.kntlim) go to 150
  130 loc=nodplc(loc)
      go to 120
  140 if (kntr.eq.0) go to 175
  150 if (ititle.eq.0) write (iofile,111)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt3) asrci,(ai1(i),i=1,kntr)
      kntr=0
      if ((loc.ne.0).and.(nodplc(loc+13).eq.0)) go to 130
c
c  nonlinear voltage controlled voltage sources
c
  175 if (jelcnt(6).eq.0) go to 186
      ititle=0
  176 format(1h0,/,'0**** voltage-controlled voltage sources')
      loc=locate(6)
      kntr=0
  178 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 182
      kntr=kntr+1
      locv=nodplc(loc+1)
      loct=lx0+nodplc(loc+13)
      anam(kntr)=value(locv)
      av1(kntr)=value(loct)
      ai1(kntr)=value(loct+1)
      if (kntr.ge.kntlim) go to 184
  180 loc=nodplc(loc)
      go to 178
  182 if (kntr.eq.0) go to 186
  184 if (ititle.eq.0) write (iofile,176)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) asrcv,(av1(i),i=1,kntr)
      write (iofile,afmt3) asrci,(ai1(i),i=1,kntr)
      kntr=0
      if ((loc.ne.0).and.(nodplc(loc+14).eq.0)) go to 180
c
c  nonlinear current controlled current sources
c
  186 if (jelcnt(7).eq.0) go to 196
      ititle=0
  187 format(1h0,/,'0**** current-controlled current sources')
      loc=locate(7)
      kntr=0
  188 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 192
      kntr=kntr+1
      locv=nodplc(loc+1)
      loct=lx0+nodplc(loc+12)
      anam(kntr)=value(locv)
      ai1(kntr)=value(loct)
      if (kntr.ge.kntlim) go to 194
  190 loc=nodplc(loc)
      go to 188
  192 if (kntr.eq.0) go to 196
  194 if (ititle.eq.0) write (iofile,187)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt3) asrci,(ai1(i),i=1,kntr)
      kntr=0
      if ((loc.ne.0).and.(nodplc(loc+13).eq.0)) go to 190
c
c  nonlinear current controlled voltage sources
c
  196 if (jelcnt(8).eq.0) go to 210
      ititle=0
  197 format(1h0,/,'0**** current-controlled voltage sources')
      loc=locate(8)
      kntr=0
  198 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 202
      kntr=kntr+1
      locv=nodplc(loc+1)
      loct=lx0+nodplc(loc+13)
      anam(kntr)=value(locv)
      av1(kntr)=value(loct)
      ai1(kntr)=value(loct+1)
      if (kntr.ge.kntlim) go to 204
  200 loc=nodplc(loc)
      go to 198
  202 if (kntr.eq.0) go to 210
  204 if (ititle.eq.0) write (iofile,197)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt2) asrcv,(av1(i),i=1,kntr)
      write (iofile,afmt3) asrci,(ai1(i),i=1,kntr)
      kntr=0
      if ((loc.ne.0).and.(nodplc(loc+14).eq.0)) go to 200
c
c  diodes
c
  210 if (jelcnt(11).eq.0) go to 300
      ititle=0
  211 format(1h0,/,'0**** diodes')
      loc=locate(11)
      kntr=0
  220 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 240
      kntr=kntr+1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+11)
      anam(kntr)=value(locv)
      amod(kntr)=value(locm)
      cd(kntr)=value(loct+1)
      vd(kntr)=value(lvnim1+node1)-value(lvnim1+node2)
      if (modedc.ne.1) go to 225
      req(kntr)=1.0d0/value(loct+2)
      cap(kntr)=value(loct+4)
  225 if (kntr.ge.kntlim) go to 250
  230 loc=nodplc(loc)
      go to 220
  240 if (kntr.eq.0) go to 300
  250 if (ititle.eq.0) write (iofile,211)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt4) (amod(i),i=1,kntr)
      write (iofile,afmt3) aid,(cd(i),i=1,kntr)
      write (iofile,afmt2) avd,(vd(i),i=1,kntr)
      if (modedc.ne.1) go to 260
      write (iofile,afmt3) areq,(req(i),i=1,kntr)
      write (iofile,afmt3) acap,(cap(i),i=1,kntr)
  260 kntr=0
      if ((loc.ne.0).and.(nodplc(loc+16).eq.0)) go to 230
c
c  bipolar junction transistors
c
  300 if (jelcnt(12).eq.0) go to 400
      ititle=0
  301 format(1h0,/,'0**** bipolar junction transistors')
      loc=locate(12)
      kntr=0
  320 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 340
      kntr=kntr+1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+8)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+22)
      anam(kntr)=value(locv)
      amod(kntr)=value(locm)
      cb(kntr)=type*value(loct+3)
      cc(kntr)=type*value(loct+2)
      vbe(kntr)=value(lvnim1+node2)-value(lvnim1+node3)
      vbc(kntr)=value(lvnim1+node2)-value(lvnim1+node1)
      vce(kntr)=vbe(kntr)-vbc(kntr)
      betadc(kntr)=cc(kntr)/dsign(dmax1(dabs(cb(kntr)),1.0d-20),
     1  cb(kntr))
      if (modedc.ne.1) go to 325
      rx(kntr)=0.0d0
      if(value(loct+16).ne.0.0d0) rx(kntr)=1.0d0/value(loct+16)
      ccs(kntr)=value(loct+13)
      cbx(kntr)=value(loct+15)
      rpi(kntr)=1.0d0/value(loct+4)
      gm(kntr)=value(loct+6)
      ro(kntr)=1.0d0/value(loct+7)
      cpi(kntr)=value(loct+9)
      cmu(kntr)=value(loct+11)
      betaac(kntr)=gm(kntr)*rpi(kntr)
      ft(kntr)=gm(kntr)/(twopi*dmax1(cpi(kntr)+cmu(kntr)+cbx(kntr),
     1  1.0d-20))
  325 if (kntr.ge.kntlim) go to 350
  330 loc=nodplc(loc)
      go to 320
  340 if (kntr.eq.0) go to 400
  350 if (ititle.eq.0) write (iofile,301)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt4) (amod(i),i=1,kntr)
      write (iofile,afmt3) aib,(cb(i),i=1,kntr)
      write (iofile,afmt3) aic,(cc(i),i=1,kntr)
      write (iofile,afmt2) avbe,(vbe(i),i=1,kntr)
      write (iofile,afmt2) avbc,(vbc(i),i=1,kntr)
      write (iofile,afmt2) avce,(vce(i),i=1,kntr)
      write (iofile,afmt2) abetad,(betadc(i),i=1,kntr)
      if (modedc.ne.1) go to 360
      write (iofile,afmt3) agm,(gm(i),i=1,kntr)
      write (iofile,afmt3) arpi,(rpi(i),i=1,kntr)
      write (iofile,afmt3) arx,(rx(i),i=1,kntr)
      write (iofile,afmt3) aro,(ro(i),i=1,kntr)
      write (iofile,afmt3) acpi,(cpi(i),i=1,kntr)
      write (iofile,afmt3) acmu,(cmu(i),i=1,kntr)
      write (iofile,afmt3) acbx,(cbx(i),i=1,kntr)
      write (iofile,afmt3) accs,(ccs(i),i=1,kntr)
      write (iofile,afmt2) abetaa,(betaac(i),i=1,kntr)
      write (iofile,afmt3) aft,(ft(i),i=1,kntr)
  360 kntr=0
      if ((loc.ne.0).and.(nodplc(loc+36).eq.0)) go to 330
c
c  jfets
c
  400 if (jelcnt(13).eq.0) go to 500
      ititle=0
  401 format(1h0,/,'0**** jfets')
      loc=locate(13)
      kntr=0
  420 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 440
      kntr=kntr+1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+7)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+19)
      anam(kntr)=value(locv)
      amod(kntr)=value(locm)
      cd(kntr)=type*(value(loct+3)-value(loct+4))
      vgs(kntr)=value(lvnim1+node2)-value(lvnim1+node3)
      vds(kntr)=value(lvnim1+node1)-value(lvnim1+node3)
      if (modedc.ne.1) go to 425
      gm(kntr)=value(loct+5)
      gds(kntr)=value(loct+6)
      cgs(kntr)=value(loct+9)
      cgd(kntr)=value(loct+11)
  425 if (kntr.ge.kntlim) go to 450
  430 loc=nodplc(loc)
      go to 420
  440 if (kntr.eq.0) go to 500
  450 if (ititle.eq.0) write (iofile,401)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt4) (amod(i),i=1,kntr)
      write (iofile,afmt3) aid,(cd(i),i=1,kntr)
      write (iofile,afmt2) avgs,(vgs(i),i=1,kntr)
      write (iofile,afmt2) avds,(vds(i),i=1,kntr)
      if (modedc.ne.1) go to 460
      write (iofile,afmt3) agm,(gm(i),i=1,kntr)
      write (iofile,afmt3) agds,(gds(i),i=1,kntr)
      write (iofile,afmt3) acgs,(cgs(i),i=1,kntr)
      write (iofile,afmt3) acgd,(cgd(i),i=1,kntr)
  460 kntr=0
      if ((loc.ne.0).and.(nodplc(loc+25).eq.0)) go to 430
c
c  mosfets
c
  500 if (jelcnt(14).eq.0) go to 600
      ititle=0
  501 format(1h0,/,'0**** mosfets')
      loc=locate(14)
      kntr=0
  520 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 540
      kntr=kntr+1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      locm=nodplc(loc+8)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+26)
      anam(kntr)=value(locv)
      amod(kntr)=value(locm)
      cd(kntr)=type*value(loct+4)
      vgs(kntr)=value(lvnim1+node2)-value(lvnim1+node3)
      vds(kntr)=value(lvnim1+node1)-value(lvnim1+node3)
      vbs(kntr)=value(lvnim1+node4)-value(lvnim1+node3)
      if (modedc.ne.1) go to 525
      xl=value(locv+1)-2.0d0*value(locm+28)
      xw=value(locv+2)
      covlgs=value(locm+13)*xw
      covlgd=value(locm+14)*xw
      covlgb=value(locm+15)*xl
      xqco=value(locm+35)
      devmod=value(locv+8)
      vdsat(kntr)=value(locv+10)
      vth(kntr)=value(locv+9)
      gm(kntr)=value(loct+7)
      gds(kntr)=value(loct+8)
      gmb(kntr)=value(loct+9)
      if(devmod.gt.0.0d0) go to 521
      vth(kntr)=value(locv+9)
  521 cbd(kntr)=value(loct+24)
      cbs(kntr)=value(loct+26)
      cgsov(kntr)=covlgs
      cgdov(kntr)=covlgd
      cgbov(kntr)=covlgb
      if (xqco.gt.0.5d0) go to 522
      cggb(kntr)=value(loct+18)
      cgdb(kntr)=value(loct+19)
      cgsb(kntr)=value(loct+20)
      cbgb(kntr)=value(loct+21)
      cbdb(kntr)=value(loct+22)
      cbsb(kntr)=value(loct+23)
      go to 525
  522 cgs(kntr)=value(loct+12)
      cgd(kntr)=value(loct+14)
      cgb(kntr)=value(loct+16)
  525 if (kntr.ge.kntlim) go to 550
  530 loc=nodplc(loc)
      go to 520
  540 if (kntr.eq.0) go to 600
  550 if (ititle.eq.0) write (iofile,501)
      ititle=1
      write (iofile,afmt1) (anam(i),i=1,kntr)
      write (iofile,afmt4) (amod(i),i=1,kntr)
      if(type.eq.0.0d0) go to 555
      write (iofile,afmt3) aid,(cd(i),i=1,kntr)
      write (iofile,afmt2) avgs,(vgs(i),i=1,kntr)
      write (iofile,afmt2) avds,(vds(i),i=1,kntr)
      write (iofile,afmt2) avbs,(vbs(i),i=1,kntr)
      if (modedc.ne.1) go to 560
      write (iofile,afmt2) avth,(vth(i),i=1,kntr)
      write (iofile,afmt2) avdsat,(vdsat(i),i=1,kntr)
      write (iofile,afmt3) agm,(gm(i),i=1,kntr)
      write (iofile,afmt3) agds,(gds(i),i=1,kntr)
      write (iofile,afmt3) agmb,(gmb(i),i=1,kntr)
      write (iofile,afmt3) acbd,(cbd(i),i=1,kntr)
      write (iofile,afmt3) acbs,(cbs(i),i=1,kntr)
      write (iofile,afmt3) acgsov,(cgsov(i),i=1,kntr)
      write (iofile,afmt3) acgdov,(cgdov(i),i=1,kntr)
      write (iofile,afmt3) acgbov,(cgbov(i),i=1,kntr)
      if (xqco.gt.0.5d0) go to 552
      write (iofile,551)
  551 format(' derivatives of gate (dqgdvx) and bulk (dqbdvx) charges')
      write (iofile,afmt3) acggb,(cggb(i),i=1,kntr)
      write (iofile,afmt3) acgdb,(cgdb(i),i=1,kntr)
      write (iofile,afmt3) acgsb,(cgsb(i),i=1,kntr)
      write (iofile,afmt3) acbgb,(cbgb(i),i=1,kntr)
      write (iofile,afmt3) acbdb,(cbdb(i),i=1,kntr)
      write (iofile,afmt3) acbsb,(cbsb(i),i=1,kntr)
      go to 560
  552 write (iofile,afmt3) acgs,(cgs(i),i=1,kntr)
      write (iofile,afmt3) acgd,(cgd(i),i=1,kntr)
      write (iofile,afmt3) acgb,(cgb(i),i=1,kntr)
      go to 560
  555 write (iofile,afmt3) aid,(cd(i),i=1,kntr)
      write (iofile,afmt3) aig,(cg(i),i=1,kntr)
      write (iofile,afmt2) avgs,(vgs(i),i=1,kntr)
      write (iofile,afmt2) avds,(vds(i),i=1,kntr)
      write (iofile,afmt2) avbs,(vbs(i),i=1,kntr)
      if (modedc.ne.1) go to 560
      write (iofile,afmt3) agm,(gm(i),i=1,kntr)
      write (iofile,afmt3) agds,(gds(i),i=1,kntr)
      write (iofile,afmt3) acgs,(cgs(i),i=1,kntr)
      write (iofile,afmt3) acgd,(cgd(i),i=1,kntr)
      write (iofile,afmt3) acgb,(cgb(i),i=1,kntr)
      write (iofile,afmt3) acds,(cds(i),i=1,kntr)
  560 kntr=0
      if ((loc.ne.0).and.(nodplc(loc+33).eq.0)) go to 530
c
c  operating point analyses
c
  600 if (modedc.ne.1) go to 700
      if (kinel.eq.0) go to 610
      call sstf
  610 if (nsens.eq.0) go to 700
      call sencal
c
c  finished
c
  700 if (modedc.eq.2) go to 710
      if (jacflg.ne.0) go to 705
      call clrmem(lvnim1)
      call clrmem(lx0)
  705 call clrmem(lvn)
      call clrmem(lvntmp)
      if (memptr(macins)) call clrmem(macins)
  710 call second(t2)
      rstats(5)=rstats(5)+t2-t1
      return
      end
