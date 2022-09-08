      subroutine disto(loco)
      implicit double precision (a-h,o-z)
c
c     this routine performs the small-signal distortion analysis.
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
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
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
      complex difvn1,difvn2,difvn3,difvi1,difvi2,difvi3,dsgo2,dsgm2,
     1   dsgmu2,dsgpi2,dscb1,dscb1r,dscje1,dscjc1,disto1,disto2,disto3,
     2   dsgmo2,dgm2o3,dgmo23,bew,cew,bcw,be2w,ce2w,bc2w,bew2,cew2,
     3   bcw2,bew12,cew12,bcw12,dscdb1,dscdj1,dsg2,cvabe,cvabc,cvace,
     4   cvout,cvdist
      dimension distit(4)
      dimension vdo(2,12)
      complex cvdo(12)
      real vdo
      equivalence (cvdo(1),vdo(1,1))
      data distit / 8hdistorti, 8hon analy, 8hsis     , 8h        /
c
c
      icvw1=ld1
      icv2w1=icvw1+nstop
      icvw2=icv2w1+nstop
      icvw12=icvw2+nstop
      icvadj=icvw12+nstop
      iprnt=0
      if (icalc.ge.2) go to 10
      idnp=nodplc(idist+2)
      idnn=nodplc(idist+3)
      locv=nodplc(idist+1)
      rload=1.0d0/value(locv+1)
      kntr=1
   10 if (idprt.eq.0) go to 30
      if (kntr.gt.icalc) go to 30
      iprnt=1
      kntr=kntr+idprt
      call title(0,lwidth,1,distit)
   30 freq1=dble(real(cvalue(loco+1)))
      freq2=skw2*freq1
      call copy16(cvalue(lcvn+1),cvalue(icvw1+1),nstop)
      cvout=cvalue(icvw1+idnp)-cvalue(icvw1+idnn)
      call magphs(cvout,omag,ophase)
c
c  begin the distortion analysis
c
      do 1000 kdisto=1,7
      cvdist=cmplx(0.0e0,0.0e0)
      go to (1000,110,120,130,140,160,170),kdisto
  110 freqd=2.0d0*freq1
      arg=dsqrt(2.0d0*rload*refprl)/(omag*omag)
      if (iprnt.eq.0) go to 200
      write (iofile,111) freq1,freqd,omag,ophase
  111 format (///5x,'2nd harmonic distortion',30x,'freq1 = ',1pd9.2,
     1   '  hz'//5x,'distortion frequency  ',d9.2,'  hz',16x,
     2   'mag ',d9.3,3x,'phs ',0pf7.2)
      go to 200
  120 freqd=3.0d0*freq1
      arg=2.0d0*rload*refprl/(omag*omag*omag)
      if (iprnt.eq.0) go to 200
      write (iofile,121) freq1,freqd,omag,ophase
  121 format (1h1,4x,'3rd harmonic distortion',30x,'freq1 = ',1pd9.2,
     1   '  hz'//5x,'distortion frequency  ',d9.2,'  hz',16x,
     2   'mag ',d9.3,3x,'phs ',0pf7.2)
      go to 200
  130 freqd=freq2
      go to 200
  140 freqd=freq1-freq2
      arg=dsqrt(2.0d0*rload*refprl)*spw2/(omag*omag)
      if (iprnt.eq.0) go to 200
      write (iofile,151) freq1,freq2,freqd,omag,ophase,ow2mag,ow2phs
  151 format (1h1,4x,'2nd order intermodulation difference component',
     1   7x,'freq1 = ',1pd9.2,'  hz',15x,'freq2 = ',d9.2,'  hz'//
     2   5x,'distortion frequency  ',d9.2,'  hz',16x,'mag ',
     3   d9.3,3x,'phs ',0pf7.2,9x,'mag ',1pd9.3,3x,'phs ',0pf7.2)
      go to 200
  160 freqd=freq1+freq2
      arg=dsqrt(2.0d0*rload*refprl)*spw2/(omag*omag)
      if (iprnt.eq.0) go to 200
      write (iofile,161) freq1,freq2,freqd,omag,ophase,ow2mag,ow2phs
  161 format (1h1,4x,'2nd order intermodulation sum component',
     1   14x,'freq1 = ',1pd9.2,'  hz',15x,'freq2 = ',d9.2,'  hz'//
     2   5x,'distortion frequency  ',d9.2,'  hz',16x,'mag ',
     3   d9.3,3x,'phs ',0pf7.2,9x,'mag ',1pd9.3,3x,'phs ',0pf7.2)
      go to 200
  170 freqd=2.0d0*freq1-freq2
      arg=2.0d0*rload*refprl*spw2/(omag*omag*omag)
      if (iprnt.eq.0) go to 200
      write (iofile,171) freq1,freq2,freqd,omag,ophase,ow2mag,ow2phs
  171 format (1h1,4x,'3rd order intermodulation difference component',
     1   7x,'freq1 = ',1pd9.2,'  hz',15x,'freq2 = ',d9.2,'  hz'//
     2   5x,'distortion frequency  ',d9.2,'  hz',16x,'mag ',
     3   d9.3,3x,'phs ',0pf7.2,9x,'mag ',1pd9.3,3x,'phs ',0pf7.2)
c
c  load and decompose y matrix
c
  200 omega=twopi*freqd
      igoof=0
      call acload
      call acdcmp
      if (igoof.eq.0) go to 220
      write (iofile,211) igoof,freqd
  211 format('0warning:  underflow ',i4,' time(s) in distortion analysis
     1 at freq = ',1pd9.3,' hz')
      igoof=0
  220 if (kdisto.eq.4) go to 710
c
c  obtain adjoint solution
c
      call zero8(value(lvn+1),nstop)
      call zero8(value(imvn+1),nstop)
      value(lvn+idnp)=-1.0d0
      value(lvn+idnn)=+1.0d0
      call acasol
      call copy16(cvalue(lcvn+1),cvalue(icvadj+1),nstop)
      call zero8(value(lvn+1),nstop)
      call zero8(value(imvn+1),nstop)
c
c  bjts
c
      if (jelcnt(12).eq.0) go to 500
      ititle=0
  301 format (////1x,'bjt distortion components'//1x,'name',11x,'gm',
     1   8x,'gpi',7x,'go',8x,'gmu',6x,'gmo2',7x,'cb',8x,'cbr',7x,'cje',
     2   7x,'cjc',6x,'total')
  311 format (////1x,'bjt distortion components'//1x,'name',11x,'gm',
     1   8x,'gpi',7x,'go',8x,'gmu',6x,'gmo2',7x,'cb',8x,'cbr',7x,'cje',
     2   7x,'cjc',6x,'gm203',5x,'gmo23',5x,'total')
  320 loc=locate(12)
  330 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 500
      locv=nodplc(loc+1)
      loct=lx0+nodplc(loc+22)
      locd=ld0+nodplc(loc+23)
      node1=nodplc(loc+5)
      node2=nodplc(loc+6)
      node3=nodplc(loc+7)
      cje1=value(locd)
      cje2=value(locd+1)
      cjc1=value(locd+2)
      cjc2=value(locd+3)
      go2=value(locd+4)
      gmo2=value(locd+5)
      gm2=value(locd+6)
      gmu2=value(locd+7)
      gpi2=value(locd+8)
      cb1=value(locd+11)
      cb1r=value(locd+12)
      go3=value(locd+13)
      gmo23=value(locd+14)
      gm2o3=value(locd+15)
      gm3=value(locd+16)
      gmu3=value(locd+17)
      gpi3=value(locd+18)
      cb2=value(locd+19)
      cb2r=value(locd+20)
      bew=cvalue(icvw1+node2)-cvalue(icvw1+node3)
      cew=cvalue(icvw1+node1)-cvalue(icvw1+node3)
      bcw=cvalue(icvw1+node2)-cvalue(icvw1+node1)
      if (kdisto.eq.2) go to 370
      be2w=cvalue(icv2w1+node2)-cvalue(icv2w1+node3)
      ce2w=cvalue(icv2w1+node1)-cvalue(icv2w1+node3)
      bc2w=cvalue(icv2w1+node2)-cvalue(icv2w1+node1)
      if (kdisto.eq.3) go to 380
      bew2=cvalue(icvw2+node2)-cvalue(icvw2+node3)
      cew2=cvalue(icvw2+node1)-cvalue(icvw2+node3)
      bcw2=cvalue(icvw2+node2)-cvalue(icvw2+node1)
      if (kdisto.eq.5) go to 390
      if (kdisto.eq.6) go to 400
      bew12=cvalue(icvw12+node2)-cvalue(icvw12+node3)
      cew12=cvalue(icvw12+node1)-cvalue(icvw12+node3)
      bcw12=cvalue(icvw12+node2)-cvalue(icvw12+node1)
      go to 410
c
c  calculate hd2 current generators
c
  370 difvn1=0.5d0*cew*cew
      difvn2=0.5d0*bew*bew
      difvn3=0.5d0*bcw*bcw
      dsgmo2=gmo2*0.5d0*bew*cew
      go to 420
c
c  calculate hd3 current generators
c
  380 difvi1=0.50d0*cew*ce2w
      difvn1=0.25d0*cew*cew*cew
      difvi2=0.50d0*bew*be2w
      difvn2=0.25d0*bew*bew*bew
      difvi3=0.50d0*bcw*bc2w
      difvn3=0.25d0*bcw*bcw*bcw
      dsgmo2=gmo2*(bew*ce2w+be2w*cew)*0.5d0
      go to 430
c
c  calculate im2d current generators
c
  390 difvn1=cew*conjg(cew2)
      difvn2=bew*conjg(bew2)
      difvn3=bcw*conjg(bcw2)
      dsgmo2=gmo2*0.5d0*(bew*conjg(cew2)+cew*conjg(bew2))
      go to 420
c
c  calculate im2s current generators
c
  400 difvn1=cew*cew2
      difvn2=bew*bew2
      difvn3=bcw*bcw2
      dsgmo2=gmo2*0.5d0*(bew*cew2+bew2*cew)
      go to 420
c
c  calculate im3 current generators
c
  410 difvi1=0.5d0*(ce2w*conjg(cew2)+cew*cew12)
      difvi2=0.5d0*(be2w*conjg(bew2)+bew*bew12)
      difvi3=0.5d0*(bc2w*conjg(bcw2)+bcw*bcw12)
      difvn1=cew*cew*conjg(cew2)*0.75d0
      difvn2=bew*bew*conjg(bew2)*0.75d0
      difvn3=bcw*bcw*conjg(bcw2)*0.75d0
      dsgmo2=gmo2*0.5d0*(conjg(bew2)*ce2w+bew*cew12+conjg(cew2)*be2w+
     1   cew*bew12)
      go to 430
c
  420 dsgo2=go2*difvn1
      dsgm2=gm2*difvn2
      dsgmu2=gmu2*difvn3
      dsgpi2=gpi2*difvn2
      dscb1=0.5d0*cb1*omega*cmplx(-aimag(difvn2),real(difvn2))
      dscb1r=0.5d0*cb1r*omega*cmplx(-aimag(difvn3),real(difvn3))
      dscje1=0.5d0*cje1*omega*cmplx(-aimag(difvn2),real(difvn2))
      dscjc1=0.5d0*cjc1*omega*cmplx(-aimag(difvn3),real(difvn3))
      go to 440
c
  430 dsgo2=2.0d0*go2*difvi1+go3*difvn1
      dsgm2=2.0d0*gm2*difvi2+gm3*difvn2
      dsgmu2=2.0d0*gmu2*difvi3+gmu3*difvn3
      dsgpi2=2.0d0*gpi2*difvi2+gpi3*difvn2
      dscb1=omega*(cb1*difvi2+cb2*difvn2/3.0d0)
      dscb1=cmplx(-aimag(dscb1),real(dscb1))
      dscb1r=omega*(cb1r*difvi3+cb2r*difvn3/3.0d0)
      dscb1r=cmplx(-aimag(dscb1r),real(dscb1r))
      dscje1=omega*(cje1*difvi2+cje2*difvn2/3.0d0)
      dscje1=cmplx(-aimag(dscje1),real(dscje1))
      dscjc1=omega*(cjc1*difvi3+cjc2*difvn3/3.0d0)
      dscjc1=cmplx(-aimag(dscjc1),real(dscjc1))
c
c  determine contribution of each distortion source
c
  440 cvabe=cvalue(icvadj+node2)-cvalue(icvadj+node3)
      cvabc=cvalue(icvadj+node2)-cvalue(icvadj+node1)
      cvace=cvalue(icvadj+node1)-cvalue(icvadj+node3)
      disto1=dsgm2+dsgo2+dsgmo2
      disto2=dsgpi2+dscb1+dscje1
      disto3=dsgmu2+dscb1r+dscjc1
      cvdo(1)=dsgm2*cvace*arg
      cvdo(2)=dsgpi2*cvabe*arg
      cvdo(3)=dsgo2*cvace*arg
      cvdo(4)=dsgmu2*cvabc*arg
      cvdo(5)=dsgmo2*cvace*arg
      cvdo(6)=dscb1*cvabe*arg
      cvdo(7)=dscb1r*cvabc*arg
      cvdo(8)=dscje1*cvabe*arg
      cvdo(9)=dscjc1*cvabc*arg
      if (kdisto.eq.3) go to 450
      if (kdisto.eq.7) go to 460
      cvdo(10)=cvdo(1)+cvdo(2)+cvdo(3)+cvdo(4)+cvdo(5)+cvdo(6)+cvdo(7)+
     1   cvdo(8)+cvdo(9)
      cvdist=cvdist+cvdo(10)
      if (iprnt.eq.0) go to 480
      do 445 j=1,10
      call magphs(cvdo(j),xmag,xphs)
      cvdo(j)=cmplx(sngl(xmag),sngl(xphs))
  445 continue
      if (ititle.eq.0) write (iofile,301)
      ititle=1
      write (iofile,446) value(locv),(vdo(1,j),j=1,10)
  446 format(1h0,a8,'mag',1p12d10.3)
      write (iofile,447) (vdo(2,j),j=1,10)
  447 format(9x,'phs',12(1x,f7.2,2x))
      go to 480
  450 dgm2o3=gm2o3*cew*bew*bew*0.25d0
      dgmo23=gmo23*bew*cew*cew*0.25d0
      go to 470
  460 dgm2o3=gm2o3*(0.5d0*bew*conjg(bew2)*cew+0.25d0*bew*bew*
     1  conjg(cew2))
      dgmo23=gmo23*(0.5d0*cew*conjg(cew2)*bew+0.25d0*cew*cew*
     1  conjg(bew2))
  470 disto1=disto1+dgm2o3+dgmo23
      cvdo(10)=dgm2o3*cvace*arg
      cvdo(11)=dgmo23*cvace*arg
      cvdo(12)=cvdo(1)+cvdo(2)+cvdo(3)+cvdo(4)+cvdo(5)+cvdo(6)+cvdo(7)+
     1   cvdo(8)+cvdo(9)+cvdo(10)+cvdo(11)
      cvdist=cvdist+cvdo(12)
      if (iprnt.eq.0) go to 480
      do 475 j=1,12
      call magphs(cvdo(j),xmag,xphs)
      cvdo(j)=cmplx(sngl(xmag),sngl(xphs))
  475 continue
      if (ititle.eq.0) write (iofile,311)
      ititle=1
      write (iofile,446) value(locv),(vdo(1,j),j=1,12)
      write (iofile,447) (vdo(2,j),j=1,12)
  480 value(lvn+node1)=value(lvn+node1)
     1  -real(disto1-disto3)
      value(lvn+node2)=value(lvn+node2)
     1  -real(disto2+disto3)
      value(lvn+node3)=value(lvn+node3)
     1  +real(disto1+disto2)
      value(imvn+node1)=value(imvn+node1)
     1  -aimag(disto1-disto3)
      value(imvn+node2)=value(imvn+node2)
     1  -aimag(disto2+disto3)
      value(imvn+node3)=value(imvn+node3)
     1  +aimag(disto1+disto2)
      loc=nodplc(loc)
      go to 330
c
c   junction diodes
c
  500 if (jelcnt(11).eq.0) go to 700
      ititle=0
  501 format (////1x,'diode distortion components'//1x,'name',
     1   11x,'geq',7x,'cb',8x,'cj',7x,'total')
  510 loc=locate(11)
  520 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 700
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      loct=lx0+nodplc(loc+11)
      locd=ld0+nodplc(loc+12)
      cdj1=value(locd)
      cdj2=value(locd+1)
      cdb1=value(locd+3)
      geq2=value(locd+4)
      geq3=value(locd+5)
      cdb2=value(locd+6)
      bew=cvalue(icvw1+node3)-cvalue(icvw1+node2)
      if (kdisto.eq.2) go to 540
      be2w=cvalue(icv2w1+node3)-cvalue(icv2w1+node2)
      if (kdisto.eq.3) go to 550
      bew2=cvalue(icvw2+node3)-cvalue(icvw2+node2)
      if (kdisto.eq.5) go to 560
      if (kdisto.eq.6) go to 570
      bew12=cvalue(icvw12+node3)-cvalue(icvw12+node2)
      go to 580
c
c    calculate hd2 current generators
c
  540 difvn1=0.5d0*bew*bew
      go to 590
c
c    calculate hd3 current generators
c
  550 difvi1=0.5d0*bew*be2w
      difvn1=0.25d0*bew*bew*bew
      go to 600
c
c    calculate im2d current generators
c
  560 difvn1=bew*conjg(bew2)
      go to 590
c
c    calculate im2s current generators
c
  570 difvn1=bew*bew2
      go to 590
c
c    calculate im3 current generators
c
  580 difvi1=0.5d0*(be2w*conjg(bew2)+bew*bew12)
      difvn1=bew*bew*conjg(bew2)*0.75d0
      go to 600
  590 dsg2=geq2*difvn1
      dscdb1=0.5d0*cdb1*omega*cmplx(-aimag(difvn1),real(difvn1))
      dscdj1=0.5d0*cdj1*omega*cmplx(-aimag(difvn1),real(difvn1))
      go to 610
c
  600 dsg2=2.0d0*geq2*difvi1+geq3*difvn1
      dscdb1=omega*(cdb1*difvi1+cdb2*difvn1/3.0d0)
      dscdb1=cmplx(-aimag(dscdb1),real(dscdb1))
      dscdj1=omega*(cdj1*difvi1+cdj2*difvn1/3.0d0)
      dscdj1=cmplx(-aimag(dscdj1),real(dscdj1))
c
c  determine contribution of each distortion source
c
  610 cvabe=cvalue(icvadj+node3)-cvalue(icvadj+node2)
      disto1=dsg2+dscdb1+dscdj1
      cvdo(1)=dsg2*cvabe*arg
      cvdo(2)=dscdb1*cvabe*arg
      cvdo(3)=dscdj1*cvabe*arg
      cvdo(4)=cvdo(1)+cvdo(2)+cvdo(3)
      cvdist=cvdist+cvdo(4)
      if (iprnt.eq.0) go to 680
      do 670 j=1,4
      call magphs(cvdo(j),xmag,xphs)
      cvdo(j)=cmplx(sngl(xmag),sngl(xphs))
  670 continue
      if (ititle.eq.0) write (iofile,501)
      ititle=1
      write (iofile,446) value(locv),(vdo(1,j),j=1,4)
      write (iofile,447) (vdo(2,j),j=1,4)
  680 value(lvn+node2)=value(lvn+node2)+real(disto1)
      value(lvn+node3)=value(lvn+node3)-real(disto1)
      value(imvn+node2)=value(imvn+node2)+aimag(disto1)
      value(imvn+node3)=value(imvn+node3)-aimag(disto1)
      loc=nodplc(loc)
      go to 520
c
c  obtain total distortion solution if necessary
c
  700 go to (1000,710,790,710,710,840,860),kdisto
  710 call acsol
c
c  store solution, print and store answers
c
  760 go to (1000,770,790,800,820,840,860),kdisto
  770 call copy16(cvalue(lcvn+1),cvalue(icv2w1+1),nstop)
      call magphs(cvdist,o2mag,o2phs)
      if (iprnt.eq.0) go to 900
      o2log=20.0d0*dlog10(o2mag)
      write (iofile,781) o2mag,o2phs,o2log
  781 format (///5x,'hd2     magnitude  ',1pd10.3,5x,'phase  ',0pf7.2,
     1   5x,'=  ',f7.2,'  db')
      go to 900
  790 call magphs(cvdist,o3mag,o3phs)
      if (iprnt.eq.0) go to 900
      o3log=20.0d0*dlog10(o3mag)
      write (iofile,791) o3mag,o3phs,o3log
  791 format (///5x,'hd3     magnitude  ',1pd10.3,5x,'phase  ',0pf7.2,
     1   5x,'=  ',f7.2,'  db')
      go to 900
  800 call copy16(cvalue(lcvn+1),cvalue(icvw2+1),nstop)
      cvout=cvalue(icvw2+idnp)-cvalue(icvw2+idnn)
      call magphs(cvout,ow2mag,ow2phs)
      go to 1000
  820 call copy16(cvalue(lcvn+1),cvalue(icvw12+1),nstop)
  840 call magphs(cvdist,o12mag,o12phs)
      if (iprnt.eq.0) go to 900
      o12log=20.0d0*dlog10(o12mag)
      if (kdisto.eq.6) go to 850
      write (iofile,841) o12mag,o12phs,o12log
  841 format (///5x,'im2d    magnitude  ',1pd10.3,5x,'phase  ',0pf7.2,
     1   5x,'=  ',f7.2,'  db')
      go to 900
  850 write (iofile,851) o12mag,o12phs,o12log
  851 format (///5x,'im2s    magnitude  ',1pd10.3,5x,'phase  ',0pf7.2,
     1   5x,'=  ',f7.2,'  db')
      go to 900
  860 call magphs(cvdist,o21mag,o21phs)
      if (iprnt.eq.0) go to 900
      o21log=20.0d0*dlog10(o21mag)
      write (iofile,861) o21mag,o21phs,o21log
  861 format (///5x,'im3     magnitude  ',1pd10.3,5x,'phase  ',0pf7.2,
     1   5x,'=  ',f7.2,'  db')
      cma=dabs(4.0d0*o21mag*dcos((o21phs-ophase)/rad))
      cma=dmax1(cma,1.0d-20)
      cmp=dabs(4.0d0*o21mag*dsin((o21phs-ophase)/rad))
      cmp=dmax1(cmp,1.0d-20)
      cmalog=20.0d0*dlog10(cma)
      cmplog=20.0d0*dlog10(cmp)
      write (iofile,866)
  866 format (////5x,'approximate cross modulation components')
      write (iofile,871) cma,cmalog
  871 format (/5x,'cma     magnitude  ',1pd10.3,24x,'=  ',0pf7.2,'  db')
      write (iofile,881) cmp,cmplog
  881 format (/5x,'cmp     magnitude  ',1pd10.3,24x,'=  ',0pf7.2,'  db')
c
c  save distortion outputs
c
  900 iflag=kdisto+2
      if (iflag.ge.7) iflag=iflag-1
      loc=locate(45)
  910 if (loc.eq.0) go to 1000
      if (nodplc(loc+5).ne.iflag) go to 920
      iseq=nodplc(loc+4)
      cvalue(loco+iseq)=cvdist
  920 loc=nodplc(loc)
      go to 910
 1000 continue
c
c  finished
c
 2000 return
      end
