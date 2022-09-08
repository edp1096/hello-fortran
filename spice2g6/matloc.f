      subroutine matloc
      implicit double precision (a-h,o-z)
c
c     this routine stores the locations of the various matrix terms to
c which the different circuit elements contribute.
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
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c  resistors
c
      loc=locate(1)
  690 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 700
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      nodplc(loc+4)=indxx(node1,node2)
      nodplc(loc+5)=indxx(node2,node1)
      nodplc(loc+6)=indxx(node1,node1)
      nodplc(loc+7)=indxx(node2,node2)
      loc=nodplc(loc)
      go to 690
c
c  capacitors
c
  700 loc=locate(2)
  710 if ((loc.eq.0).or.(nodplc(loc+12).ne.0)) go to 720
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      nodplc(loc+5)=indxx(node1,node2)
      nodplc(loc+6)=indxx(node2,node1)
      nodplc(loc+10)=indxx(node1,node1)
      nodplc(loc+11)=indxx(node2,node2)
      loc=nodplc(loc)
      go to 710
c
c  inductors
c
  720 loc=locate(3)
  730 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 740
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ibr=nodplc(loc+5)
      nodplc(loc+6)=indxx(node1,ibr)
      nodplc(loc+7)=indxx(node2,ibr)
      nodplc(loc+8)=indxx(ibr,node1)
      nodplc(loc+9)=indxx(ibr,node2)
      nodplc(loc+13)=indxx(ibr,ibr)
      loc=nodplc(loc)
      go to 730
c
c  mutual inductances
c
  740 loc=locate(4)
  750 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 760
      nl1=nodplc(loc+2)
      nl2=nodplc(loc+3)
      ibr1=nodplc(nl1+5)
      ibr2=nodplc(nl2+5)
      nodplc(loc+4)=indxx(ibr1,ibr2)
      nodplc(loc+5)=indxx(ibr2,ibr1)
      loc=nodplc(loc)
      go to 750
c
c  nonlinear voltage controlled current sources
c
  760 loc=locate(5)
  762 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 764
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      lnod=nodplc(loc+6)
      lmat=nodplc(loc+7)
      do 763 i=1,ndim
      node3=nodplc(lnod+1)
      node4=nodplc(lnod+2)
      lnod=lnod+2
      nodplc(lmat+1)=indxx(node1,node3)
      nodplc(lmat+2)=indxx(node1,node4)
      nodplc(lmat+3)=indxx(node2,node3)
      nodplc(lmat+4)=indxx(node2,node4)
      lmat=lmat+4
  763 continue
      loc=nodplc(loc)
      go to 762
c
c  nonlinear voltage controlled voltage sources
c
  764 loc=locate(6)
  766 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 768
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      ibr=nodplc(loc+6)
      lnod=nodplc(loc+7)
      lmat=nodplc(loc+8)
      nodplc(lmat+1)=indxx(node1,ibr)
      nodplc(lmat+2)=indxx(node2,ibr)
      nodplc(lmat+3)=indxx(ibr,node1)
      nodplc(lmat+4)=indxx(ibr,node2)
      lmat=lmat+4
      do 767 i=1,ndim
      node3=nodplc(lnod+1)
      node4=nodplc(lnod+2)
      lnod=lnod+2
      nodplc(lmat+1)=indxx(ibr,node3)
      nodplc(lmat+2)=indxx(ibr,node4)
      lmat=lmat+2
  767 continue
      loc=nodplc(loc)
      go to 766
c
c  nonlinear current controlled current sources
c
  768 loc=locate(7)
  770 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 772
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      locvs=nodplc(loc+6)
      lmat=nodplc(loc+7)
      do 771 i=1,ndim
      locvst=nodplc(locvs+i)
      ibr=nodplc(locvst+6)
      nodplc(lmat+1)=indxx(node1,ibr)
      nodplc(lmat+2)=indxx(node2,ibr)
      lmat=lmat+2
  771 continue
      loc=nodplc(loc)
      go to 770
c
c  nonlinear current controlled voltage sources
c
  772 loc=locate(8)
  774 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 780
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      ndim=nodplc(loc+4)
      ibr=nodplc(loc+6)
      locvs=nodplc(loc+7)
      lmat=nodplc(loc+8)
      nodplc(lmat+1)=indxx(node1,ibr)
      nodplc(lmat+2)=indxx(node2,ibr)
      nodplc(lmat+3)=indxx(ibr,node1)
      nodplc(lmat+4)=indxx(ibr,node2)
      lmat=lmat+4
      do 775 i=1,ndim
      locvst=nodplc(locvs+i)
      kbr=nodplc(locvst+6)
      nodplc(lmat+i)=indxx(ibr,kbr)
  775 continue
      loc=nodplc(loc)
      go to 774
c
c  voltage sources
c
  780 loc=locate(9)
  790 if ((loc.eq.0).or.(nodplc(loc+11).ne.0)) go to 800
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      iptr=nodplc(loc+6)
      nodplc(loc+7)=indxx(node1,iptr)
      nodplc(loc+8)=indxx(node2,iptr)
      nodplc(loc+9)=indxx(iptr,node1)
      nodplc(loc+10)=indxx(iptr,node2)
      loc=nodplc(loc)
      go to 790
c
c  diodes
c
  800 loc=locate(11)
  810 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 820
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      nodplc(loc+7)=indxx(node1,node3)
      nodplc(loc+8)=indxx(node2,node3)
      nodplc(loc+9)=indxx(node3,node1)
      nodplc(loc+10)=indxx(node3,node2)
      nodplc(loc+13)=indxx(node1,node1)
      nodplc(loc+14)=indxx(node2,node2)
      nodplc(loc+15)=indxx(node3,node3)
      loc=nodplc(loc)
      go to 810
c
c  transistors
c
  820 loc=locate(12)
  830 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 840
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      node7=nodplc(loc+30)
      nodplc(loc+10)=indxx(node1,node4)
      nodplc(loc+11)=indxx(node2,node5)
      nodplc(loc+12)=indxx(node3,node6)
      nodplc(loc+13)=indxx(node4,node1)
      nodplc(loc+14)=indxx(node4,node5)
      nodplc(loc+15)=indxx(node4,node6)
      nodplc(loc+16)=indxx(node5,node2)
      nodplc(loc+17)=indxx(node5,node4)
      nodplc(loc+18)=indxx(node5,node6)
      nodplc(loc+19)=indxx(node6,node3)
      nodplc(loc+20)=indxx(node6,node4)
      nodplc(loc+21)=indxx(node6,node5)
      nodplc(loc+24)=indxx(node1,node1)
      nodplc(loc+25)=indxx(node2,node2)
      nodplc(loc+26)=indxx(node3,node3)
      nodplc(loc+27)=indxx(node4,node4)
      nodplc(loc+28)=indxx(node5,node5)
      nodplc(loc+29)=indxx(node6,node6)
      nodplc(loc+31)=indxx(node7,node7)
      nodplc(loc+32)=indxx(node4,node7)
      nodplc(loc+33)=indxx(node7,node4)
      nodplc(loc+34)=indxx(node2,node4)
      nodplc(loc+35)=indxx(node4,node2)
      loc=nodplc(loc)
      go to 830
c
c  jfets
c
  840 loc=locate(13)
  850 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 860
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      nodplc(loc+9)=indxx(node1,node4)
      nodplc(loc+10)=indxx(node2,node4)
      nodplc(loc+11)=indxx(node2,node5)
      nodplc(loc+12)=indxx(node3,node5)
      nodplc(loc+13)=indxx(node4,node1)
      nodplc(loc+14)=indxx(node4,node2)
      nodplc(loc+15)=indxx(node4,node5)
      nodplc(loc+16)=indxx(node5,node2)
      nodplc(loc+17)=indxx(node5,node3)
      nodplc(loc+18)=indxx(node5,node4)
      nodplc(loc+20)=indxx(node1,node1)
      nodplc(loc+21)=indxx(node2,node2)
      nodplc(loc+22)=indxx(node3,node3)
      nodplc(loc+23)=indxx(node4,node4)
      nodplc(loc+24)=indxx(node5,node5)
      loc=nodplc(loc)
      go to 850
c
c  mosfets
c
  860 loc=locate(14)
  870 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 900
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      node5=nodplc(loc+6)
      node6=nodplc(loc+7)
      nodplc(loc+10)=indxx(node1,node5)
      nodplc(loc+11)=indxx(node2,node4)
      nodplc(loc+12)=indxx(node2,node5)
      nodplc(loc+13)=indxx(node2,node6)
      nodplc(loc+14)=indxx(node3,node6)
      nodplc(loc+15)=indxx(node4,node2)
      nodplc(loc+16)=indxx(node4,node5)
      nodplc(loc+17)=indxx(node4,node6)
      nodplc(loc+18)=indxx(node5,node1)
      nodplc(loc+19)=indxx(node5,node2)
      nodplc(loc+20)=indxx(node5,node4)
      nodplc(loc+21)=indxx(node5,node6)
      nodplc(loc+22)=indxx(node6,node2)
      nodplc(loc+23)=indxx(node6,node3)
      nodplc(loc+24)=indxx(node6,node4)
      nodplc(loc+25)=indxx(node6,node5)
      nodplc(loc+27)=indxx(node1,node1)
      nodplc(loc+28)=indxx(node2,node2)
      nodplc(loc+29)=indxx(node3,node3)
      nodplc(loc+30)=indxx(node4,node4)
      nodplc(loc+31)=indxx(node5,node5)
      nodplc(loc+32)=indxx(node6,node6)
      loc=nodplc(loc)
      go to 870
c
c  transmission lines
c
  900 loc=locate(17)
  910 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 1000
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      ni1=nodplc(loc+6)
      ni2=nodplc(loc+7)
      ibr1=nodplc(loc+8)
      ibr2=nodplc(loc+9)
      nodplc(loc+10)=indxx(node1,node1)
      nodplc(loc+11)=indxx(node1,ni1)
      nodplc(loc+12)=indxx(node2,ibr1)
      nodplc(loc+13)=indxx(node3,node3)
      nodplc(loc+14)=indxx(node4,ibr2)
      nodplc(loc+15)=indxx(ni1,node1)
      nodplc(loc+16)=indxx(ni1,ni1)
      nodplc(loc+17)=indxx(ni1,ibr1)
      nodplc(loc+18)=indxx(ni2,ni2)
      nodplc(loc+19)=indxx(ni2,ibr2)
      nodplc(loc+20)=indxx(ibr1,node2)
      nodplc(loc+21)=indxx(ibr1,node3)
      nodplc(loc+22)=indxx(ibr1,node4)
      nodplc(loc+23)=indxx(ibr1,ni1)
      nodplc(loc+24)=indxx(ibr1,ibr2)
      nodplc(loc+25)=indxx(ibr2,node1)
      nodplc(loc+26)=indxx(ibr2,node2)
      nodplc(loc+27)=indxx(ibr2,node4)
      nodplc(loc+28)=indxx(ibr2,ni2)
      nodplc(loc+29)=indxx(ibr2,ibr1)
      nodplc(loc+31)=indxx(node3,ni2)
      nodplc(loc+32)=indxx(ni2,node3)
      loc=nodplc(loc)
      go to 910
c
c  finished
c
 1000 return
      end
