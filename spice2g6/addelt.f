      subroutine addelt(loce,loc,id,inodx,inodi,nnodi)
      implicit double precision (a-h,o-z)
c
c     this routine adds an element to the nominal circuit definition
c lists.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c... inodx(1), inodi(1) are arrays (see subckt)
      dimension inodx(1),inodi(1)
c
      dimension lnod(50),lval(50),nnods(50)
      data lnod /10,14,16, 8,15,16,15,16,13, 8,
     1           18,38,27,35, 8, 8,35, 5, 5, 5,
     2            5, 5, 5, 5, 0, 0, 0, 0, 0, 0,
     3           21,21,21,21,21,21,21,21,21,21,
     4            8, 8, 8, 8, 8, 0, 0, 0, 0, 0 /
      data lval / 5, 4, 4, 2, 1, 1, 1, 1, 4, 4,
     1            3, 4, 4,16, 1, 1, 9, 2, 1, 1,
     2           19,55,17,46, 0, 0, 0, 0, 0, 0,
     3            1, 1, 1, 1, 1,17,17,17,17,17,
     4            1, 1, 1, 1, 1, 0, 0, 0, 0, 0 /
      data nnods / 2, 2, 2, 0, 2, 2, 2, 2, 2, 2,
     1             2, 4, 3, 4, 4, 4, 4, 0, 1, 0,
     2             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     3             0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     4             2, 2, 2, 0, 0, 0, 0, 0, 0, 0 /
c
c  copy integer part
c
      nword=lnod(id)-3
      if (nword.le.0) go to 10
      call copy4(nodplc(loc+2),nodplc(loce+2),nword)
c
c  set nodes
c
   10 if (id.ge.21) go to 100
      if (nnods(id).eq.0) go to 100
      if (id.le.4) go to 20
      if (id.le.8) go to 40
      if (id.eq.19) go to 70
   20 jstop=nnods(id)
      do 30 j=1,jstop
      call newnod(nodplc(loc+j+1),nodplc(loce+j+1),inodx(1),
     1  inodi(1),nnodi)
   30 continue
      go to 100
   40 call newnod(nodplc(loc+2),nodplc(loce+2),inodx(1),inodi(1),nnodi)
      call newnod(nodplc(loc+3),nodplc(loce+3),inodx(1),inodi(1),nnodi)
      if (id.ge.7) go to 100
      nlocp=loc+id+1
      nssnod=2*nodplc(loc+4)
      call getm4(nodplc(loce+id+1),nssnod)
      nlocpe=loce+id+1
   50 do 60 j=1,nssnod
      locp=nodplc(nlocp)
      nodold=nodplc(locp+j)
      call newnod(nodold,nodnew,inodx(1),inodi(1),nnodi)
      locpe=nodplc(nlocpe)
      nodplc(locpe+j)=nodnew
   60 continue
      go to 100
   70 nlocp=loc+2
      call sizmem(nodplc(loc+2),nssnod)
      call getm4(nodplc(loce+2),nssnod)
      nlocpe=loce+2
      go to 50
c
c  copy real part
c
  100 if (nogo.ne.0) go to 300
      locv=nodplc(loc+1)
      locve=nodplc(loce+1)
      call copy8(value(locv),value(locve),lval(id))
c
c  treat non-node tables specially
c
  200 if (id.ge.11) go to 300
      go to (300,210,220,300,230,240,230,240,260,260), id
  210 if (nodplc(loc+4).eq.1) go to 300
      call cpytb8(loc+7,loce+7)
      go to 300
  220 if (nodplc(loc+4).eq.1) go to 300
      call cpytb8(loc+10,loce+10)
      go to 300
  230 itab=5
      go to 250
  240 itab=6
  250 if (id.le.6) go to 255
      call cpytb4(loc+itab+1,loce+itab+1)
  255 call cpytb4(loc+itab+2,loce+itab+2)
      call cpytb8(loc+itab+3,loce+itab+3)
      call cpytb8(loc+itab+4,loce+itab+4)
      call cpytb4(loc+itab+5,loce+itab+5)
      call cpytb8(loc+itab+6,loce+itab+6)
      go to 300
  260 call cpytb8(loc+5,loce+5)
c
c
  300 return
      end
