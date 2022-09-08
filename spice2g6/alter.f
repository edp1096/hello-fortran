      subroutine alter
      implicit double precision (a-h,o-z)
c
c     this routine changes the element or device parameters
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
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
      logical memptr
c
      integer xxor
      dimension lnod(50),lval(50)
      dimension chtitl(4)
      data chtitl / 8hchange f,8hollowing,8h paramet,8hers     /
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
c
      call title (0,lwidth,1,chtitl)
      do 350 id=1,24
      loc=locate(id)
   10 if (loc.eq.0) go to 350
      if (nodplc(loc+lnod(id)-2).ne.numcyc) go to 300
      locv=nodplc(loc+1)
      loc1=locate(id)
   50 if (loc1.eq.0) go to 400
      if (nodplc(loc1+lnod(id)-2).ne.0) go to 400
      locv1=nodplc(loc1+1)
      if (xxor(value(locv),value(locv1)).eq.0) go to 100
      loc1=nodplc(loc1)
      go to 50
c
c  copy changed values to the original tables
c
c  copy real part
c
  100 call copy8(value(locv),value(locv1),lval(id))
      write (iofile,110) value(locv1)
  110 format ('********      ',a8,'      ********')
c
c  treat non-node tables specially
c
  200 if (id.ge.11) go to 300
      go to (300,210,220,300,230,240,230,240,260,260), id
  210 if (nodplc(loc+4).eq.1) go to 300
      if (memptr(nodplc(loc1+7))) call clrmem(nodplc(loc1+7))
      call cpytb8(loc+7,loc1+7)
      go to 300
  220 if (nodplc(loc+4).eq.1) go to 300
      if (memptr(nodplc(loc1+10))) call clrmem(nodplc(loc1+10))
      call cpytb8(loc+10,loc1+10)
      go to 300
  230 itab=5
      go to 250
  240 itab=6
  250 if (id.le.6) go to 255
      if (memptr(nodplc(loc1+itab+1))) call clrmem(nodplc(loc1+itab+1))
      call cpytb4(loc+itab+1,loc1+itab+1)
  255 if (memptr(nodplc(loc1+itab+2))) call clrmem(nodplc(loc1+itab+2))
      call cpytb4(loc+itab+2,loc1+itab+2)
      if (memptr(nodplc(loc1+itab+3))) call clrmem(nodplc(loc1+itab+3))
      call cpytb8(loc+itab+3,loc1+itab+3)
      if (memptr(nodplc(loc1+itab+4))) call clrmem(nodplc(loc1+itab+4))
      call cpytb8(loc+itab+4,loc1+itab+4)
      if (memptr(nodplc(loc1+itab+5))) call clrmem(nodplc(loc1+itab+5))
      call cpytb4(loc+itab+5,loc1+itab+5)
      if (memptr(nodplc(loc1+itab+6))) call clrmem(nodplc(loc1+itab+6))
      call cpytb8(loc+itab+6,loc1+itab+6)
      go to 300
  260 if (memptr(nodplc(loc1+5))) call clrmem(nodplc(loc1+5))
      call cpytb8(loc+5,loc1+5)
c
  300 loc=nodplc(loc)
      go to 10
  350 continue
      write (iofile,360)
  360 format (//)
      go to 500
c
  400 write (iofile,401) value(nodplc(loc1+1))
  401 format ('0*error*:  parameter change failed',/,
     1        '0*******:  ',a8,' is not in the original circuit')
      nogo=1
c
  500 return
      end
