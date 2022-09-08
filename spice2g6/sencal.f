      subroutine sencal
      implicit double precision (a-h,o-z)
c
c     this routine computes the dc sensitivities of circuit elements
c with respect to user specified outputs.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension string(5),sentit(4)
      data alsrs,alsis,alsn,alsrb,alsrc,alsre / 2hrs,2his,1hn,2hrb,2hrc,
     1   2hre /
      data alsbf,alsise,alsbr,alsisc,alsne,alsnc,alsik,alsikr,alsva,alsvb
     1   / 2hbf,3hise,2hbr,3hisc,2hne,2hnc,3hikf,3hikr,3hvaf,3hvar/
      data alsjs /2hjs/
      data sentit / 8hdc sensi, 8htivity a, 8hnalysis , 8h         /
      data ablnk / 1h  /
c
c
      if (kinel.ne.0) go to 8
    4 call dcdcmp
c
c
    8 do 1000 n=1,nsens
c
c  prepare adjoint excitation vector
c
      call zero8(value(lvn+1),nstop)
      locs=nodplc(isens+n)
      ioutyp=nodplc(locs+5)
      if (ioutyp.ne.0) go to 10
c...  voltage output
      ivolts=1
      noposo=nodplc(locs+2)
      nonego=nodplc(locs+3)
      value(lvn+noposo)=-1.0d0
      value(lvn+nonego)=+1.0d0
      go to 20
c...  current output (through voltage source)
   10 iptro=nodplc(locs+2)
      ivolts=0
      iptro=nodplc(iptro+6)
      value(lvn+iptro)=-1.0d0
c
c  obtain adjoint solution by doing forward/backward substitution on
c  the transpose of the y matrix
c
   20 call asol
      value(lvn+1)=0.0d0
c
c  real solution in lvnim1;  adjoint solution in lvn ...
c
      call title(0,lwidth,1,sentit)
      ipos=1
      call outnam(locs,1,string,ipos)
      call move(string,ipos,ablnk,1,7)
      jstop=(ipos+6)/8
      write (iofile,36) (string(j),j=1,jstop)
   36 format('0dc sensitivities of output ',5a8)
      if(ivolts.ne.0) write (iofile,41)
      if(ivolts.eq.0) write(iofile,42)
   41 format(1h0,8x,'element',9x,'element',7x,'element',7x,'normalized'/
     1   10x,'name',12x,'value',6x,'sensitivity    sensitivity'/35x,
     2   ' (volts/unit) (volts/percent)'/)
   42 format(1h0,8x,'element',9x,'element',7x,'element',7x,'normalized'/
     1   10x,'name',12x,'value',6x,'sensitivity    sensitivity'/35x,
     2   '  (amps/unit)  (amps/percent)'/)
c
c  resistors
c
      loc=locate(1)
  100 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 110
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      val=1.0d0/value(locv+1)
      sens=-(value(lvnim1+node1)-value(lvnim1+node2))*
     1      (value(lvn   +node1)-value(lvn   +node2))/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) value(locv),val,sens,sensn
  101 format(10x,a8,4x,1pd10.3,5x,d10.3,5x,d10.3)
  105 loc=nodplc(loc)
      go to 100
c
c  voltage sources
c
  110 loc=locate(9)
  140 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 150
      locv=nodplc(loc+1)
      val=value(locv+1)
      iptrv=nodplc(loc+6)
      sens=-value(lvn+iptrv)
      sensn=val*sens/100.0d0
      write (iofile,101) value(locv),val,sens,sensn
  145 loc=nodplc(loc)
      go to 140
c
c  current sources
c
  150 loc=locate(10)
  160 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 170
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      val=value(locv+1)
      sens=value(lvn+node1)-value(lvn+node2)
      sensn=val*sens/100.0d0
      write (iofile,101) value(locv),val,sens,sensn
  165 loc=nodplc(loc)
      go to 160
c
c  diodes
c
  170 loc=locate(11)
  180 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 210
      locv=nodplc(loc+1)
      write (iofile,181) value(locv)
  181 format(1x,a8)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      area=value(locv+1)
c
c  series resistance (rs)
c
      val=value(locm+2)*area
      if (val.ne.0.0d0) go to 190
      write (iofile,186) alsrs
  186 format(10x,a8,5x,2h0.,13x,2h0.,13x,2h0.)
      go to 200
  190 val=1.0d0/val
      sens=-(value(lvnim1+node1)-value(lvnim1+node3))*
     1      (value(lvn   +node1)-value(lvn   +node3))/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsrs,val,sens,sensn
c
c  intrinsic parameters
c
  200 csat=value(locm+1)*area
      xn=value(locm+3)
      vbe=value(lvnim1+node3)-value(lvnim1+node2)
      vte=xn*vt
      evbe=dexp(vbe/vte)
      vabe=value(lvn+node3)-value(lvn+node2)
c
c  saturation current (is)
c
      sens=vabe*(evbe-1.0d0)
      sensn=csat*sens/100.0d0
      write (iofile,101) alsis,csat,sens,sensn
c
c  ideality factor (n)
c
      sens=-vabe*(csat/xn)*(vbe/vte)*evbe
      if (dabs(sens).lt.1.0d-30) sens=0.0d0
      sensn=xn*sens/100.0d0
      write (iofile,101) alsn,xn,sens,sensn
  205 loc=nodplc(loc)
      go to 180
c
c  bipolar junction transistors
c
  210 loc=locate(12)
  220 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 1000
      locv=nodplc(loc+1)
      write (iofile,181) value(locv)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      locm=nodplc(loc+8)
      type=nodplc(locm+2)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+22)
      area=value(locv+1)
c
c  base resistance (rb)
c
      val=value(loct+16)
      if (val.ne.0.0d0) go to 230
      write (iofile,186) alsrb
      go to 240
  230 val=1.0d0/val
      sens=-(value(lvnim1+node2)-value(lvnim1+node5))*
     1      (value(lvn   +node2)-value(lvn   +node5))/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsrb,val,sens,sensn
c
c  collector resistance (rc)
c
  240 val=value(locm+20)*area
      if (val.ne.0.0d0) go to 250
      write (iofile,186) alsrc
      go to 260
  250 val=1.0d0/val
      sens=-(value(lvnim1+node1)-value(lvnim1+node4))*
     1      (value(lvn   +node1)-value(lvn   +node4))/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsrc,val,sens,sensn
c
c  emitter resistance (re)
c
  260 val=value(locm+19)*area
      if (val.ne.0.0d0) go to 270
      write (iofile,186) alsre
      go to 280
  270 val=1.0d0/val
      sens=-(value(lvnim1+node3)-value(lvnim1+node6))*
     1      (value(lvn   +node3)-value(lvn   +node6))/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsre,val,sens,sensn
c
c  intrinsic parameters
c
  280 bf=value(locm+2)
      br=value(locm+8)
      csat=value(locm+1)*area
      ova=value(locm+4)
      ovb=value(locm+10)
      oik=value(locm+5)/area
      ise=value(locm+6)*area
      xne=value(locm+7)
      vte=xne*vt
      oikr=value(locm+11)/area
      isc=value(locm+12)*area
      xnc=value(locm+13)
      vtc=xnc*vt
      vbe=type*(value(lvnim1+node5)-value(lvnim1+node6))
      vbc=type*(value(lvnim1+node5)-value(lvnim1+node4))
      vabe=type*(value(lvn+node5)-value(lvn+node6))
      vabc=type*(value(lvn+node5)-value(lvn+node4))
      vace=vabe-vabc
      if (vbe.le.-vt) go to 320
      evbe=dexp(vbe/vt/value(locm+3))
      cbe=csat*(evbe-1.0d0)
      gbe=csat*evbe/vt/value(locm+3)
      if (ise.ne.0.0d0) go to 310
      cben=0.0d0
      gben=0.0d0
      go to 350
  310 evben=dexp(vbe/vte)
      cben=ise     *(evben-1.0d0)
      gben=ise     *evben/vte
      go to 350
  320 gbe=-csat/vbe
      cbe=gbe*vbe
      gben=-ise/vbe
      cben=gben*vbe
  350 if (vbc.le.-vt) go to 370
      evbc=dexp(vbc/vt/value(locm+9))
      cbc=csat*(evbc-1.0d0)
      gbc=csat*evbc/vt/value(locm+9)
      if (isc.ne.0.0d0) go to 360
      cbcn=0.0d0
      gbcn=0.0d0
      go to 400
  360 evbcn=dexp(vbc/vtc)
      cbcn=isc     *(evbcn-1.0d0)
      gbcn=isc     *evbcn/vtc
      go to 400
  370 gbc=-csat/vbc
      cbc=gbc*vbc
      gbcn=-isc/vbc
      cbcn=gbcn*vbc
  400 q1=1.0d0/(1.0d0-ova*vbc-ovb*vbe)
      q2=oik*cbe+oikr*cbc
      sqarg=dsqrt(1.0d0+4.0d0*q2)
      qb=q1*(1.0d0+sqarg)/2.0d0
      dqb=(cbe-cbc)/(qb*qb)
      sqarg=dsqrt(1.0d0+4.0d0*q2)
      dq1=dqb*(1.0d0+sqarg)/2.0d0
      dq2=q1*dqb/sqarg
c
c  compute sensitivities
c
c...  bf
      sens=-vabe*cbe/bf/bf
      sensn=bf*sens/100.0d0
      write (iofile,101) alsbf,bf,sens,sensn
c...  ise
      if (ise.ne.0.0d0) go to 430
      write (iofile,186) alsise
      go to 440
  430 sens=vabe*cben/ise
      sensn=ise*sens/100.0d0
      write (iofile,101) alsise,ise,sens,sensn
c...  br
  440 sens=-vabc*cbc/br/br
      sensn=br*sens/100.0d0
      write (iofile,101) alsbr,br,sens,sensn
c...  isc
      if (isc.ne.0.0d0) go to 450
      write (iofile,186) alsisc
      go to 460
  450 sens=vabc*cbcn/isc
      sensn=isc*sens/100.0d0
      write (iofile,101) alsisc,isc,sens,sensn
c...  is
  460 sens=(vabe*(cbe/bf)+vabc*(cbc/br)
     1   +vace*(dqb*qb-dq2*q2))/csat
      sensn=csat*sens/100.0d0
      write (iofile,101) alsjs,csat,sens,sensn
c...  ne
      sens=-vabe*gben*vbe/xne
      sensn=xne*sens/100.0d0
      write (iofile,101) alsne,xne,sens,sensn
c...  nc
      sens=-vabc*gbcn*vbc/xnc
      sensn=xnc*sens/100.0d0
      write (iofile,101) alsnc,xnc,sens,sensn
c...  ik
      if (oik.ne.0.0d0) go to 470
      write (iofile,186) alsik
      go to 480
  470 val=1.0d0/oik
      sens=vace*dq2*cbe/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsik,val,sens,sensn
c...  ikr
  480 if (oikr.ne.0.0d0) go to 490
      write (iofile,186) alsikr
      go to 500
  490 val=1.0d0/oikr
      sens=vace*dq2*cbc/(val*val)
      sensn=val*sens/100.0d0
      write (iofile,101) alsikr,val,sens,sensn
c...  va
  500 if (ova.ne.0.0d0) go to 510
      write (iofile,186) alsva
      go to 520
  510 va=1.0d0/ova
      sens=vace*q1*q1*dq1*vbc/(va*va)
      sensn=va*sens/100.0d0
      write (iofile,101) alsva,va,sens,sensn
c...  vb
  520 if (ovb.ne.0.0d0) go to 530
      write (iofile,186) alsvb
      go to 540
  530 vb=1.0d0/ovb
      sens=vace*q1*q1*dq1*vbe/(vb*vb)
      sensn=vb*sens/100.0d0
      write (iofile,101) alsvb,vb,sens,sensn
c
c
  540 loc=nodplc(loc)
      go to 220
c
c  finished
c
 1000 continue
      return
      end
