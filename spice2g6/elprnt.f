      subroutine elprnt
      implicit double precision (a-h,o-z)
c
c     this routine prints a circuit element summary.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension itab(25),astyp(6)
      dimension eltitl(4)
      data eltitl / 8hcircuit , 8helement , 8hsummary , 8h         /
      data astyp / 1h , 5hpulse, 3hsin, 3hexp, 3hpwl, 4hsffm /
      data ablnk,aoff /1h ,3hoff/
c
c  print listing of elements
c
      call title(0,lwidth,1,eltitl)
c
c  print resistors
c
      if (jelcnt(1).eq.0) go to 50
      ititle=0
   21 format(//'0**** resistors'/'0     name        nodes     value
     1  tc1        tc2'//)
      loc=locate(1)
   30 if ((loc.eq.0).or.(nodplc(loc+8).ne.0)) go to 50
      if (ititle.eq.0) write (iofile,21)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      write (iofile,31) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),value(locv+2),value(locv+3),value(locv+4)
   31 format(6x,a8,2i5,1p3d11.2)
   40 loc=nodplc(loc)
      go to 30
c
c  print capacitors and inductors
c
   50 if ((jelcnt(2)+jelcnt(3)).eq.0) go to 80
      ititle=0
   51 format(//'0**** capacitors and inductors'/'0     name        nodes
     1    in cond     value'//)
      do 70 id=2,3
      loc=locate(id)
   60 if (loc.eq.0) go to 70
      if ((id.eq.2).and.(nodplc(loc+12).ne.0)) go to 70
      if ((id.eq.3).and.(nodplc(loc+14).ne.0)) go to 70
      if (ititle.eq.0) write (iofile,51)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      if (nodplc(loc+4).ne.1) go to 62
      write (iofile,31) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),value(locv+2),value(locv+1)
      go to 65
   62 ltab=7
      if (id.eq.3) ltab=10
      call sizmem(nodplc(loc+ltab),nparam)
      ispot=nodplc(loc+ltab)+1
      write (iofile,63) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),value(locv+2)
   63 format(6x,a8,2i5,1pd11.2,'   variable')
   65 loc=nodplc(loc)
      go to 60
   70 continue
c
c  print mutual inductors
c
   80 if (jelcnt(4).eq.0) go to 100
      ititle=0
   81 format(//'0**** mutual inductors'/'0     name        coupled induc
     1tors   value'//)
      loc=locate(4)
   90 if ((loc.eq.0).or.(nodplc(loc+6).ne.0)) go to 110
      if (ititle.eq.0) write (iofile,81)
      ititle=1
      locv=nodplc(loc+1)
      nl1=nodplc(loc+2)
      nl1=nodplc(nl1+1)
      nl2=nodplc(loc+3)
      nl2=nodplc(nl2+1)
      write (iofile,91) value(locv),value(nl1),value(nl2),value(locv+1)
   91 format(6x,a8,4x,a8,2x,a8,1pd10.2)
   95 loc=nodplc(loc)
      go to 90
c
c  print nonlinear voltage controlled sources
c
  100 if (jelcnt(5).eq.0) go to 120
      ititle=0
  101 format(//'0**** voltage-controlled current sources'/'0     name
     1     +    -   dimension   function')
      loc=locate(5)
  110 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 120
      if (ititle.eq.0) write (iofile,101)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      write (iofile,111) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(loc+4)
  111 format(6x,a8,2i5,i8,9x,'poly')
  115 loc=nodplc(loc)
      go to 110
c
c  nonlinear voltage controlled voltage sources
c
  120 if (jelcnt(6).eq.0) go to 140
      ititle=0
  121 format(//'0**** voltage-controlled voltage sources'/'0     name
     1     +    -   dimension   function')
      loc=locate(6)
  130 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 140
      if (ititle.eq.0) write (iofile,121)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      write (iofile,111) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(loc+4)
  135 loc=nodplc(loc)
      go to 130
c
c  nonlinear current controlled current sources
c
  140 if (jelcnt(7).eq.0) go to 160
      ititle=0
  141 format(//'0**** current-controlled current sources'/'0     name
     1     +    -   dimension   function')
      loc=locate(7)
  150 if ((loc.eq.0).or.(nodplc(loc+13).ne.0)) go to 160
      if (ititle.eq.0) write (iofile,141)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      write (iofile,111) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(loc+4)
  155 loc=nodplc(loc)
      go to 150
c
c  nonlinear current controlled voltage sources
c
  160 if (jelcnt(8).eq.0) go to 170
      ititle=0
  161 format(//'0**** current-controlled voltage sources'/'0     name
     1     +    -   dimension   function')
      loc=locate(8)
  165 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 170
      if (ititle.eq.0) write (iofile,161)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      write (iofile,111) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(loc+4)
  167 loc=nodplc(loc)
      go to 165
c
c  print independent sources
c
  170 if ((jelcnt(9)+jelcnt(10)).eq.0) go to 250
      ititle=0
  171 format(//'0**** independent sources'/'0     name        nodes   dc
     1 value   ac value   ac phase   transient'//)
      do 245 id=9,10
      loc=locate(id)
  180 if (loc.eq.0) go to 245
      if ((id.eq.9).and.(nodplc(loc+11).ne.0)) go to 245
      if ((id.eq.10).and.(nodplc(loc+6).ne.0)) go to 245
      if (ititle.eq.0) write (iofile,171)
      ititle=1
      locv=nodplc(loc+1)
      locp=nodplc(loc+5)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      itype=nodplc(loc+4)+1
      anam=astyp(itype)
      write (iofile,181) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),value(locv+1),value(locv+2),
     2   value(locv+3),anam
  181 format(6x,a8,2i5,1p3d11.2,2x,a8)
      if (jtrflg.eq.0) go to 240
      jstart=locp+1
      go to (240,190,200,210,220,230), itype
  190 jstop=locp+7
      write (iofile,191) (value(j),j=jstart,jstop)
  191 format(1h0,42x,'initial value',1pd11.2,/,
     1           43x,'pulsed value.',  d11.2,/,
     2           43x,'delay time...',  d11.2,/,
     3           43x,'risetime.....',  d11.2,/,
     4           43x,'falltime.....',  d11.2,/,
     5           43x,'width........',  d11.2,/,
     6           43x,'period.......',  d11.2,/)
      go to 240
  200 jstop=locp+5
      write (iofile,201) (value(j),j=jstart,jstop)
  201 format(1h0,42x,'offset.......',1pd11.2,/,
     1           43x,'amplitude....',  d11.2,/,
     2           43x,'frequency....',  d11.2,/,
     3           43x,'delay........',  d11.2,/,
     4           43x,'theta........',  d11.2,/)
      go to 240
  210 jstop=locp+6
      write (iofile,211) (value(j),j=jstart,jstop)
  211 format(1h0,42x,'initial value',1pd11.2,/,
     1           43x,'pulsed value.',  d11.2,/,
     2           43x,'rise delay...',  d11.2,/,
     3           43x,'rise tau.....',  d11.2,/,
     4           43x,'fall delay...',  d11.2,/,
     5           43x,'fall tau.....',  d11.2,/)
      go to 240
  220 call sizmem(nodplc(loc+5),jstop)
      jstop=locp+jstop
      write (iofile,221) (value(j),j=jstart,jstop)
  221 format(1h0,49x,'time       value'//,(46x,1p2d11.2))
      write (iofile,226)
  226 format(1x)
      go to 240
  230 jstop=locp+5
      write (iofile,231) (value(j),j=jstart,jstop)
  231 format(1h0,42x,'offset.......',1pd11.2,/,
     1           43x,'amplitude....',  d11.2,/,
     2           43x,'carrier freq.',  d11.2,/,
     3           43x,'modn index...',  d11.2,/,
     4           43x,'signal freq..',  d11.2,/)
  240 loc=nodplc(loc)
      go to 180
  245 continue
c
c  print transmission lines
c
  250 if (jelcnt(17).eq.0) go to 260
      ititle=0
  251 format(//'0**** transmission lines'/'0     name             nodes
     1           z0         td'//)
      loc=locate(17)
  253 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 260
      if (ititle.eq.0) write (iofile,251)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      write (iofile,256) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(junode+node3),
     2   nodplc(junode+node4),value(locv+1),value(locv+2)
  256 format(6x,a8,4i5,1p2d11.2)
  258 loc=nodplc(loc)
      go to 253
c
c  print diodes
c
  260 if (jelcnt(11).eq.0) go to 290
      ititle=0
  261 format(//'0**** diodes'/'0     name        +    -  model       are
     1a'//)
      loc=locate(11)
  270 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 290
      if (ititle.eq.0) write (iofile,261)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      locm=nodplc(loc+5)
      locm=nodplc(locm+1)
      aic=ablnk
      if (nodplc(loc+6).eq.1) aic=aoff
      write (iofile,271) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),value(locm),value(locv+1),aic
  271 format(6x,a8,2i5,2x,a8,f8.3,2x,a8)
  280 loc=nodplc(loc)
      go to 270
c
c  print transistors
c
  290 if (jelcnt(12).eq.0) go to 320
      ititle=0
  291 format(//'0**** bipolar junction transistors'/'0     name        c
     1    b    e    s  model       area'//)
      loc=locate(12)
  300 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 320
      if (ititle.eq.0) write (iofile,291)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      aic=ablnk
      if (nodplc(loc+9).eq.1) aic=aoff
      write (iofile,301) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(junode+node3),nodplc(junode+node4),
     2   value(locm),value(locv+1),aic
  301 format(6x,a8,4i5,2x,a8,f8.3,2x,a8)
  310 loc=nodplc(loc)
      go to 300
c
c  print jfets
c
  320 if (jelcnt(13).eq.0) go to 350
      ititle=0
  321 format(//'0**** jfets'/'0     name        d    g    s  model
     1 area'//)
      loc=locate(13)
  330 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 350
      if (ititle.eq.0) write (iofile,321)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      locm=nodplc(loc+7)
      locm=nodplc(locm+1)
      aic=ablnk
      if (nodplc(loc+8).eq.1) aic=aoff
      write (iofile,331) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(junode+node3),
     2   value(locm),value(locv+1),aic
  331 format(6x,a8,3i5,2x,a8,f8.3,2x,a8)
  340 loc=nodplc(loc)
      go to 330
c
c  print mosfets
c
  350 if (jelcnt(14).eq.0) go to 400
      ititle=0
  351 format(//'0**** mosfets',/,'0name',6x,'d   g   s   b  model',6x,
     1      'w       ad       pd      rds'/
     2  37x,'l       as       ps      rss',//)
      loc=locate(14)
  360 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 400
      if (ititle.eq.0) write (iofile,351)
      ititle=1
      locv=nodplc(loc+1)
      node1=nodplc(loc+2)
      node2=nodplc(loc+3)
      node3=nodplc(loc+4)
      node4=nodplc(loc+5)
      locm=nodplc(loc+8)
      locm=nodplc(locm+1)
      aic=ablnk
      if (nodplc(loc+9).eq.1) aic=aoff
      write (iofile,361) value(locv),nodplc(junode+node1),
     1   nodplc(junode+node2),nodplc(junode+node3),
     2   nodplc(junode+node4),value(locm),value(locv+2),
     3   value(locv+3),value(locv+11),value(locv+13),
     4   value(locv+1),value(locv+4),value(locv+12),value(locv+14),aic
  361 format(1x,a8,4i4,1x,a8,1p4d8.1,/34x,1p4d8.1,1x,a8)
  370 loc=nodplc(loc)
      go to 360
c
c  subcircuit calls
c
  400 if (jelcnt(19).eq.0) go to 500
      ititle=0
  401 format(//'0**** subcircuit calls'/'0     name     subcircuit   ext
     1ernal nodes'//)
      loc=locate(19)
  410 if (loc.eq.0) go to 500
      if (ititle.eq.0) write (iofile,401)
      ititle=1
      locv=nodplc(loc+1)
      locn=nodplc(loc+2)
      call sizmem(nodplc(loc+2),nnodx)
      locs=nodplc(loc+3)
      locsv=nodplc(locs+1)
      jstart=1
      ndprln=(lwidth-28)/5
  412 jstop=min0(nnodx,jstart+ndprln-1)
      do 414 j=jstart,jstop
      node=nodplc(locn+j)
      itab(j-jstart+1)=nodplc(junode+node)
  414 continue
      if (jstart.eq.1)
     1   write (iofile,416) value(locv),value(locsv),(itab(j),j=1,jstop)
  416 format(6x,a8,2x,a8,4x,20i5)
      if (jstart.ne.1)
     1   write (iofile,418) (itab(j-jstart+1),j=jstart,jstop)
  418 format(28x,20i5)
      jstart=jstop+1
      if (jstart.le.nnodx) go to 412
      if (nnodx.le.ndprln) go to 420
      write (iofile,226)
  420 loc=nodplc(loc)
      go to 410
c
c  finished
c
  500 return
      end
