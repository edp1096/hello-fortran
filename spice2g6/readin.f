c spice version 2g.6  sccsid=readin.ma 3/22/83
      subroutine readin
      implicit double precision (a-h,o-z)
c
c
c     this routine drives the input processing of spice.  element cards
c and device models are handled by this routine.
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
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
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
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c  control card identifiers
c
      dimension aide(20),nnods(20),ntnods(20),l2nod(24)
      dimension numic(4)
      dimension aidm(7),ipolar(7),modid(7),ipar(5),ampar(115)
      dimension titinp(4)
      dimension aidc(22)
      data titinp / 8hinput li, 8hsting   , 8h        , 8h         /
      data naidc / 22 /
      data aidc / 8hac      , 8hdc      , 8hdistorti, 8hend     ,
     1            8hends    , 8hfourier , 8hmodel   , 8hnoise   ,
     2            8hop      , 8hoptions , 8hplot    , 8hprint   ,
     3            8hsubckt  , 8hsensitiv, 8htransien, 8htf      ,
     4            8htemperat, 8hwidth   , 8hnodeset , 8hic      ,
     5            8h:debug: , 8halter   /
c
c  element card identifiers, keywords, and information
c
      data aide / 1hr,1hc,1hl,1hk,1hg,1he,1hf,1hh,1hv,1hi,1hd,1hq,1hj,
     1   1hm,1hs,1hy,1ht,0.0d0,1hx,0.0d0 /
      data alsac,alspu,alsex,alssi /2hac,2hpu,2hex,2hsi/
      data alsoff,alsdc,alspw / 3hoff,2hdc,3hpw  /
      data alsz0,alszo,alsnl,alsf,alstd / 2hz0,2hzo,2hnl,1hf,2htd /
      data alsl,alsw,alsas,alsad,alspd,alsps,alsrds,alsrss,alsxqc
     1   /1hl,1hw,2has,2had,2hpd,2hps,3hnrd,3hnrs,3hxqc/
      data alszx /2hzx/
      data alssf / 4hsf   /
      data apoly, aic, area / 4hpoly, 2hic, 4harea /
      data alstc / 2htc /
      data numic / 1, 2, 2, 3 /
      data ablnk, aper / 1h , 1h. /
      data nnods / 2,2,2,0,2,2,2,2,2,2,2,3,3,4,4,4,4,0,0,0 /
      data ntnods / 2,2,2,0,2,2,2,2,2,2,3,6,5,6,4,4,4,0,0,0 /
      data l2nod / 8,12,14, 6,13,14,13,14,11, 6,
     1            16,36,25,33, 6, 6,33, 0, 3, 3,
     2             3, 3, 3, 3 /
c
c  model card keywords
c
      data aidm /1hd,3hnpn,3hpnp,3hnjf,3hpjf,4hnmos,4hpmos/
      data ipolar /0,1,-1,1,-1,1,-1/
      data modid /1,2,2,3,3,4,4/
      data ipar / 0, 14, 60, 72, 114/
      data ampar /
     1   6his    ,6hrs    ,6hn     ,6htt    ,6hcjo   ,6hvj    ,6hm     ,
     2   6heg    ,6hxti   ,6hkf    ,6haf    ,6hfc    ,6hbv    ,6hibv   ,
     1   6his    ,6hbf    ,6hnf    ,6hvaf   ,6hikf   ,6hise   ,6hne    ,
     2   6hbr    ,6hnr    ,6hvar   ,6hikr   ,6hisc   ,6hnc    ,6h0     ,
     3   6h0     ,6hrb    ,6hirb   ,6hrbm   ,6hre    ,6hrc    ,6hcje   ,
     4   6hvje   ,6hmje   ,6htf    ,6hxtf   ,6hvtf   ,6hitf   ,6hptf   ,
     5   6hcjc   ,6hvjc   ,6hmjc   ,6hxcjc  ,6htr    ,6h0     ,6h0     ,
     6   6h0     ,6h0     ,6hcjs   ,6hvjs   ,6hmjs   ,6hxtb   ,6heg    ,
     7   6hxti   ,6hkf    ,6haf    ,6hfc    ,
     1   6hvto   ,6hbeta  ,6hlambda,6hrd    ,6hrs    ,6hcgs   ,6hcgd   ,
     2   6hpb    ,6his    ,6hkf    ,6haf    ,6hfc    ,
     1   6hlevel ,6hvto   ,6hkp    ,6hgamma ,6hphi   ,6hlambda,6hrd    ,
     2   6hrs    ,6hcbd   ,6hcbs   ,6his    ,6hpb    ,6hcgso  ,6hcgdo  ,
     3   6hcgbo  ,6hrsh   ,6hcj    ,6hmj    ,6hcjsw  ,6hmjsw  ,6hjs    ,
     4   6htox   ,6hnsub  ,6hnss   ,6hnfs   ,6htpg   ,6hxj    ,6hld    ,
     5   6huo    ,6hucrit ,6huexp  ,6hutra  ,6hvmax  ,6hneff  ,6hxqc   ,
     6   6hkf    ,6haf    ,6hfc    ,6hdelta ,6htheta ,6heta   ,6hkappa ,
     7   0.0d0   /
c
c  initialize variables
c
      call second(t1)
      call zero4(idebug,20)
      call getlin
      if (keof.ne.0) go to 6000
      call copy8(afield,atitle,10)
      call getm4(ielmnt,0)
      call getm8(itemps,1)
      value(itemps+1)=27.0d0
      itemno=1
      nopage=0
      call title(-1,72,1,titinp)
      iwidth=80
      do 5 i=1,8
      achar=ablnk
      call move(achar,1,atitle(10),i,1)
      if(achar.eq.ablnk) go to 8
    5 continue
      write(iofile,6)
    6 format('0warning:  input line-width set to 72 columns because',/
     11x,'possible sequencing appears in cols 73-80')
      iwidth=72
    8 do 10 i=1,15
      afield(i)=ablnk
   10 continue
      call copy8(afield,oldlin,15)
      call getm4(isbckt,0)
      nsbckt=0
      call getm8(iunsat,0)
      nunsat=0
      numalt=0
      numcyc=0
      lwidth=132
      iprnta=0
      iprntl=0
      iprntm=1
      iprntn=0
      iprnto=0
      gmin=1.0d-12
      pivtol=1.0d-13
      pivrel=1.0d-3
      reltol=0.001d0
      abstol=1.0d-12
      vntol=1.0d-6
      trtol=7.0d0
      chgtol=1.0d-14
      defl=1.0d-4
      defw=1.0d-4
      defad=0.0d0
      defas=0.0d0
      numdgt=4
      numtem=1
      itl1=100
      itl2=50
      itl3=4
      itl4=10
      itl5=5000
      itl6=0
      limtim=2
      limpts=201
      lvlcod=1
      lvltim=2
      method=1
      xmu=0.5d0
      maxord=2
      nosolv=0
      icvflg=0
      itcelm(2)=0
      idist=0
      idprt=0
      inoise=0
      jacflg=0
      jtrflg=0
      call getm4(ifour,0)
      nfour=0
      call getm4(nsnod,0)
      call getm8(nsval,0)
      call getm4(icnod,0)
      call getm8(icval,0)
      kinel=0
      kovar=0
      kssop=0
      nosprt=0
      nsens=0
      call getm4(isens,0)
      numnod=0
      ncnods=0
      nunods=0
      call zero4(locate,50)
      call zero4(jelcnt,50)
      insize=50
      call getm8(ifield,insize)
      call getm4(icode,insize)
      call getm8(idelim,insize)
      call getm4(icolum,insize)
      go to 50
c
c  error entry
c
   40 nogo=1
c
c  read and decode next card in input deck
c
   50 igoof=0
      call card
      if (keof.ne.0) go to 5000
      if (igoof.ne.0) go to 40
      if (nodplc(icode+1).eq.0) go to 95
      anam=value(ifield+1)
      call move(anam,2,ablnk,1,7)
      if (anam.ne.aper) go to 70
      call move(anam,1,value(ifield+1),2,7)
      call keysrc(aidc,naidc,anam,id)
      if (id.le.0) go to 90
      if (id.eq.4) go to 5000
      if (id.eq.5) go to 800
      if (id.eq.7) go to 500
      if (id.eq.13) go to 700
      if (id.eq.22) numalt=numalt+1
      if (nsbckt.ge.1) go to 85
      if (id.ne.22) call runcon(id)
      if (igoof.ne.0) go to 40
      go to 50
   70 id=0
   80 id=id+1
      if (id.gt.20) go to 90
      if (anam.eq.aide(id)) go to 100
      go to 80
   85 write (iofile,86)
   86 format('0warning:  above line not allowed within subcircuit -- ',
     1   'ignored'/)
      go to 50
   90 write (iofile,91) value(ifield+1)
   91 format('0*error*:  unknown data card:  ',a8/)
      go to 40
   95 write (iofile,96)
   96 format('0*error*:  unrecognizable data card'/)
      go to 40
c
c  element and device cards
c
  100 call find(value(ifield+1),id,loc,1)
      locv=nodplc(loc+1)
      if (id.eq.4) go to 140
      if (id.eq.19) go to 900
      istop=nnods(id)+1
      if (nodplc(loc+l2nod(id)).ne.0) go to 113
      do 110 i=2,istop
      if (nodplc(icode+i).ne.0) go to 410
      if (value(ifield+i).lt.0.0d0) go to 400
  110 nodplc(loc+i)=value(ifield+i)
      go to 115
  113 do 114 i=2,istop
      nodplc(loc+i)=0
  114 continue
  115 go to (120,130,130,140,150,150,180,180,200,200,300,300,300,300,
     1   390,390,350,390,390,390), id
c
c  resistor
c
  120 if (nodplc(icode+4).ne.0) go to 420
      if (value(ifield+4).eq.0.0d0) go to 480
      value(locv+2)=value(ifield+4)
      ifld=4
  122 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,122,124
  124 anam=value(ifield+ifld)
      if (anam.ne.alstc) go to 460
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,126,124
  126 value(locv+3)=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,128,124
  128 value(locv+4)=value(ifield+ifld)
      go to 50
c
c  capacitor or inductor
c
  130 iknt=0
      ltab=7
      if (id.eq.3) ltab=10
      if (nodplc(icode+4)) 420,131,132
  131 if (value(ifield+4).le.0.0d0) go to 420
      value(locv+1)=value(ifield+4)
      nodplc(loc+4)=1
      ifld=5
      if (nodplc(icode+ifld)) 50,420,139
  132 call getm8(nodplc(loc+ltab),0)
      anam=value(ifield+4)
      if (anam.ne.apoly) go to 450
      ifld=4
  134 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,136,138
  136 call extmem(nodplc(loc+ltab),1)
      iknt=iknt+1
      ispot=nodplc(loc+ltab)+iknt
      value(ispot)=value(ifield+ifld)
      go to 134
  138 if (iknt.eq.0) go to 420
  139 anam=value(ifield+ifld)
      if (anam.ne.aic) go to 460
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 50
      value(locv+2)=value(ifield+ifld)
      go to 50
c
c  mutual inductance
c
  140 if (nodplc(icode+2).ne.1) go to 430
      anam=value(ifield+2)
      call move(anam,2,ablnk,1,7)
      if (anam.ne.aide(3)) go to 430
      call extnam(value(ifield+2),nodplc(loc+2))
      if (nodplc(icode+3).ne.1) go to 430
      anam=value(ifield+3)
      call move(anam,2,ablnk,1,7)
      if (anam.ne.aide(3)) go to 430
      call extnam(value(ifield+3),nodplc(loc+3))
      if (nodplc(icode+4).ne.0) go to 420
      xk=value(ifield+4)
      if (xk.le.0.0d0) go to 420
      if (xk.le.1.0d0) go to 145
      xk=1.0d0
      write (iofile,141)
  141 format('0warning:  coefficient of coupling reset to 1.0d0'/)
  145 value(locv+1)=xk
      go to 50
c
c  voltage controlled (nonlinear) sources
c
  150 ndim=1
      ifld=3
      if (nodplc(icode+4)) 410,156,152
  152 anam=value(ifield+4)
      if (anam.ne.apoly) go to 450
      if (nodplc(icode+5).ne.0) go to 420
      ndim=value(ifield+5)
      if (ndim.le.0) go to 420
      ifld=5
  156 nodplc(loc+4)=ndim
      ltab=id+1
      nssnod=2*ndim
      nmat=4*ndim
      if (id.eq.6) nmat=4+2*ndim
      call getm4(nodplc(loc+ltab),nssnod)
      call getm4(nodplc(loc+ltab+1),nmat)
      call getm8(nodplc(loc+ltab+2),0)
      call getm8(nodplc(loc+ltab+3),ndim)
      call getm4(nodplc(loc+ltab+4),ndim)
      call getm8(nodplc(loc+ltab+5),ndim)
      ispot=nodplc(loc+ltab+5)
      call zero8(value(ispot+1),ndim)
      lnod=nodplc(loc+ltab)
      do 158 i=1,nssnod
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 410
      if (value(ifield+ifld).lt.0.0d0) go to 400
      nodplc(lnod+i)=value(ifield+ifld)
  158 continue
  160 iknt=0
  162 ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 164
      call extmem(nodplc(loc+ltab+2),1)
      iknt=iknt+1
      ispot=nodplc(loc+ltab+2)+iknt
      value(ispot)=value(ifield+ifld)
      go to 162
  164 if (iknt.eq.0) go to 420
      if (nodplc(icode+ifld).ne.1) go to 170
      anam=value(ifield+ifld)
      if (anam.ne.aic) go to 460
      do 168 i=1,ndim
      ifld=ifld+1
      if (nodplc(icode+ifld)) 170,166,420
  166 ispot=nodplc(loc+ltab+5)+i
      value(ispot)=value(ifield+ifld)
  168 continue
  170 if (ndim.ne.1) go to 50
      if (iknt.ne.1) go to 50
      call extmem(nodplc(loc+ltab+2),1)
      ispot=nodplc(loc+ltab+2)
      value(ispot+2)=value(ispot+1)
      value(ispot+1)=0.0d0
      go to 50
c
c  current controlled (nonlinear) sources
c
  180 ndim=1
      ifld=3
      if (nodplc(icode+4).ne.1) go to 470
      anam=value(ifield+4)
      if (anam.ne.apoly) go to 182
      ifld=5
      if (nodplc(icode+5).ne.0) go to 420
      ndim=value(ifield+5)
      if (ndim.le.0) go to 420
  182 nodplc(loc+4)=ndim
      ltab=id-1
      nmat=2*ndim
      if (id.eq.8) nmat=4+ndim
      call getm4(nodplc(loc+ltab),ndim)
      call getm4(nodplc(loc+ltab+1),nmat)
      call getm8(nodplc(loc+ltab+2),0)
      call getm8(nodplc(loc+ltab+3),ndim)
      call getm4(nodplc(loc+ltab+4),ndim)
      call getm8(nodplc(loc+ltab+5),ndim)
      ispot=nodplc(loc+ltab+5)
      call zero8(value(ispot+1),ndim)
      do 184 i=1,ndim
      ifld=ifld+1
      if (nodplc(icode+ifld).ne.1) go to 470
      anam=value(ifield+ifld)
      call move(anam,2,ablnk,1,7)
      if (anam.ne.aide(9)) go to 470
      call extnam(value(ifield+ifld),loct)
      ispot=nodplc(loc+ltab)+i
      nodplc(ispot)=loct
  184 continue
      go to 160
c
c  independent sources
c
  200 ifld=3
      call getm8(nodplc(loc+5),0)
  210 ifld=ifld+1
  215 if (nodplc(icode+ifld)) 50,220,230
  220 if (ifld.gt.4) go to 210
  225 value(locv+1)=value(ifield+ifld)
      go to 210
  230 anam=value(ifield+ifld)
      if (anam.ne.alsdc) go to 235
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,225,230
  235 if (anam.ne.alsac) go to 260
      value(locv+2)=1.0d0
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,240,230
  240 value(locv+2)=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,250,230
  250 value(locv+3)=value(ifield+ifld)
      go to 210
  260 id=0
      call move(anam,3,ablnk,1,6)
      if (anam.eq.alspu) id=1
      if (anam.eq.alssi) id=2
      if (anam.eq.alsex) id=3
      if (anam.eq.alspw) id=4
      if (anam.eq.alssf) id=5
      if (id.eq.0) go to 450
      nodplc(loc+4)=id
      iknt=0
  270 ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 280
      call extmem(nodplc(loc+5),1)
      iknt=iknt+1
      ispot=nodplc(loc+5)+iknt
      value(ispot)=value(ifield+ifld)
      go to 270
  280 aval=0.0d0
      if (id.ne.4) go to 285
c...  for pwl source function, force even number of input values
      ibit=0
      if(iknt.ne.(iknt/2)*2) ibit=1
      aval=value(ispot)
      if (ibit.eq.0) go to 290
      call extmem(nodplc(loc+5),1)
      aval=value(ispot-1)
      iknt=iknt+1
      ispot=nodplc(loc+5)+iknt
      value(ispot)=aval
      go to 290
  285 if (iknt.ge.7) go to 215
  290 call extmem(nodplc(loc+5),2)
      ispot=nodplc(loc+5)+iknt
      value(ispot+1)=0.0d0
      value(ispot+2)=aval
      iknt=iknt+2
      go to 285
c
c  device cards
c
  300 value(locv+1)=1.0d0
      if (id.ne.14) go to 305
      value(locv+1)=0.0d0
      value(locv+11)=0.0d0
      value(locv+12)=0.0d0
      value(locv+13)=1.0d0
      value(locv+14)=1.0d0
      value(locv+15)=0.0d0
  305 locm=loc+ntnods(id)+2
      ifld=nnods(id)+2
c
c  temporarily (until modchk) put bjt's substrate node into nodplc(loc+5)
c
      if(id.ne.12) go to 308
      if(nodplc(icode+5).ne.0) go to 308
      ifld=6
      if (nodplc(loc+l2nod(id)).ne.0) go to 306
      nodplc(loc+5)=value(ifield+5)
      go to 308
  306 nodplc(loc+5)=0
  308 continue
c
c    reserve device internal nodes,read device geometry parameters
c    and initial conditions
c
      if (nodplc(icode+ifld).ne.1) go to 440
      call extnam(value(ifield+ifld),nodplc(locm))
  310 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,325,315
  315 anam=value(ifield+ifld)
      if (anam.ne.alsoff) go to 320
      nodplc(locm+1)=1
      go to 310
  320 if (anam.ne.area) go to 330
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,325,315
  325 if (value(ifield+ifld).le.0.0d0) go to 420
      if (id.eq.14) go to 343
      value(locv+1)=value(ifield+ifld)
      go to 310
  330 if (anam.ne.aic) go to 341
      iknt=0
      icloc=0
      if (id.eq.14) icloc=3
      maxknt=numic(id-10)
  335 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,340,315
  340 iknt=iknt+1
      if (iknt.gt.maxknt) go to 335
      value(locv+icloc+iknt+1)=value(ifield+ifld)
      go to 335
  341 if (id.ne.14) go to 460
      ispot=0
      if (anam.eq.alsl) ispot=1
      if (anam.eq.alsw) ispot=2
      if (anam.eq.alsad) ispot=3
      if (anam.eq.alszx) ispot=3
      if (anam.eq.alsas) ispot=4
      if (anam.eq.alspd) ispot=11
      if (anam.eq.alsps) ispot=12
      if (anam.eq.alsrds) ispot=13
      if (anam.eq.alsrss) ispot=14
      if (anam.eq.alsxqc) ispot=15
      if (ispot.eq.0) go to 460
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,342,315
  342 if (value(ifield+ifld).le.0.0d0) go to 420
      value(locv+ispot)=value(ifield+ifld)
      go to 310
  343 iknt=0
  344 iknt=iknt+1
      if(value(ifield+ifld).le.0.0d0) go to 420
      if(iknt.gt.15) go to 490
      if(iknt.eq.5) iknt=11
      value(locv+iknt)=value(ifield+ifld)
      ifld=ifld+1
      if(nodplc(icode+ifld)) 345,344,345
  345 if(nodplc(icode+ifld)) 50,50,315
c
c  transmission lines
c
  350 ifld=5
      xnl=0.25d0
      tfreq=0.0d0
  355 ifld=ifld+1
      if (nodplc(icode+ifld)) 378,355,360
  360 anam=value(ifield+ifld)
      if (anam.eq.aic) go to 364
      if (anam.eq.alsnl) go to 370
      if (anam.eq.alsf) go to 374
      id=0
      if (anam.eq.alsz0) id=1
      if (anam.eq.alszo) id=1
      if (anam.eq.alstd) id=2
      if (id.eq.0) go to 460
      ifld=ifld+1
      if (nodplc(icode+ifld)) 378,362,360
  362 if (value(ifield+ifld).le.0.0d0) go to 420
      value(locv+id)=value(ifield+ifld)
      go to 355
  364 iknt=0
  366 ifld=ifld+1
      if (nodplc(icode+ifld)) 378,368,360
  368 iknt=iknt+1
      if (iknt.gt.4) go to 366
      value(locv+iknt+4)=value(ifield+ifld)
      go to 366
  370 ifld=ifld+1
      if (nodplc(icode+ifld)) 378,372,360
  372 if (value(ifield+ifld).le.0.0d0) go to 420
      xnl=value(ifield+ifld)
      go to 355
  374 ifld=ifld+1
      if (nodplc(icode+ifld)) 378,376,360
  376 if (value(ifield+ifld).le.0.0d0) go to 420
      tfreq=value(ifield+ifld)
      go to 355
  378 if (value(locv+1).ne.0.0d0) go to 380
      write (iofile,379)
  379 format('0*error*:  z0 must be specified'/)
      go to 40
  380 if (value(locv+2).ne.0.0d0) go to 50
      if (tfreq.ne.0.0d0) go to 382
      write (iofile,381)
  381 format('0*error*:  either td or f must be specified'/)
      go to 40
  382 value(locv+2)=xnl/tfreq
      go to 50
c
c  elements not yet implemented
c
  390 write (iofile,391)
  391 format('0*error*:  element type not yet implemented'/)
      go to 40
c
c  element card errors
c
  400 write (iofile,401)
  401 format('0*error*:  negative node number found'/)
      go to 40
  410 write (iofile,411)
  411 format('0*error*:  node numbers are missing'/)
      go to 40
  420 write (iofile,421)
  421 format('0*error*:  value is missing or is nonpositive'/)
      go to 40
  430 write (iofile,431)
  431 format('0*error*:  mutual inductance references are missing'/)
      go to 40
  440 write (iofile,441)
  441 format('0*error*:  model name is missing'/)
      go to 40
  450 write (iofile,451) anam
  451 format('0*error*:  unknown source function:  ',a8)
      go to 40
  460 write (iofile,461) anam
  461 format('0*error*:  unknown parameter:  ',a8/)
      go to 40
  470 write (iofile,471)
  471 format('0*error*:  voltage source not found on above line'/)
      go to 40
  480 write (iofile,481)
  481 format('0*error*:  value is zero'/)
      go to 40
  490 write(iofile,491)
  491 format('0*error*:  extra numerical data on mosfet card'/)
      go to 40
c
c  model card
c
  500 if (nodplc(icode+2).ne.1) go to 650
      if (nodplc(icode+3).ne.1) go to 650
      id=0
  510 id=id+1
      if (id.gt.7) go to 660
      if (value(ifield+3).ne.aidm(id)) go to 510
      ipol=ipolar(id)
      jtype=modid(id)
      id=jtype+20
      call find(value(ifield+2),id,loc,1)
      nodplc(loc+2)=ipol
      locv=nodplc(loc+1)
  520 locm=ipar(jtype)
      nopar=ipar(jtype+1)-locm
      ifld=3
  530 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,530,560
  560 anam=value(ifield+ifld)
      if(jtype.eq.2) anam=alias(anam)
      iknt=0
  570 iknt=iknt+1
      if (iknt.gt.nopar) go to 670
      if (anam.ne.ampar(locm+iknt)) go to 570
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,580,560
  580 value(locv+iknt)=value(ifield+ifld)
      ifld=ifld+1
      if (nodplc(icode+ifld)) 50,590,560
  590 iknt=iknt+1
      if (iknt.gt.nopar) go to 530
      if (ablnk.ne.ampar(locm+iknt)) go to 530
      go to 580
c
c  model card errors
c
  650 write (iofile,651)
  651 format('0*error*:  model type is missing'/)
      go to 40
  660 write (iofile,661) value(ifield+3)
  661 format('0*error*:  unknown model type:  ',a8/)
      go to 40
  670 write (iofile,671) anam
  671 format('0*error*:  unknown model parameter:  ',a8,/)
      nogo=1
      go to 530
c
c  subcircuit definition
c
  700 if (nodplc(icode+2).ne.1) go to 780
      call find(value(ifield+2),20,loc,1)
      call extmem(isbckt,1)
      nsbckt=nsbckt+1
      nodplc(isbckt+nsbckt)=loc
      ifld=2
      if (nodplc(icode+3).ne.0) go to 790
      call getm4(nodplc(loc+2),0)
      iknt=0
  710 ifld=ifld+1
      if (nodplc(icode+ifld)) 50,720,710
  720 call extmem(nodplc(loc+2),1)
      iknt=iknt+1
      ispot=nodplc(loc+2)+iknt
      if (value(ifield+ifld).le.0.0d0) go to 770
      nodplc(ispot)=value(ifield+ifld)
      node=nodplc(ispot)
      i=iknt-1
  730 if (i.eq.0) go to 710
      ispot=ispot-1
      if (nodplc(ispot).eq.node) go to 760
      i=i-1
      go to 730
  760 write (iofile,761) node
  761 format('0*error*:  subcircuit definition duplicates node ',i5,/)
      go to 40
  770 write (iofile,771)
  771 format('0*error*:  nonpositive node number found in subcircuit ',
     1   'definition'/)
      go to 40
  780 write (iofile,781)
  781 format('0*error*:  subcircuit name missing'/)
      go to 40
  790 write (iofile,791)
  791 format('0*error*:  subcircuit nodes missing'/)
      go to 40
c
c  .ends processing
c
  800 if (nsbckt.eq.0) go to 890
      iknt=1
      if (nodplc(icode+2).le.0) go to 820
      anam=value(ifield+2)
      iknt=nsbckt
  810 loc=nodplc(isbckt+iknt)
      locv=nodplc(loc+1)
      anams=value(locv)
      if (anam.eq.anams) go to 820
      iknt=iknt-1
      if (iknt.ne.0) go to 810
      go to 880
  820 irel=nsbckt-iknt+1
      call relmem(isbckt,irel)
      nsbckt=nsbckt-irel
      go to 50
  880 write (iofile,881) anam
  881 format('0*error*:  unknown subcircuit name:  ',a8/)
      go to 40
  890 write (iofile,891)
  891 format('0warning:  no subcircuit definition known -- line ignored'
     1/)
      go to 50
c
c  subcircuit call
c
  900 call getm4(nodplc(loc+2),0)
      ifld=1
      iknt=0
  910 ifld=ifld+1
      if (nodplc(icode+ifld).ne.0) go to 920
      call extmem(nodplc(loc+2),1)
      iknt=iknt+1
      ispot=nodplc(loc+2)+iknt
      if (value(ifield+ifld).lt.0.0d0) go to 400
      nodplc(ispot)=value(ifield+ifld)
      go to 910
  920 if (iknt.eq.0) go to 410
      if (nodplc(icode+ifld).ne.1) go to 990
      call extnam(value(ifield+ifld),nodplc(loc+3))
      go to 50
  990 write (iofile,991)
  991 format('0*error*:  subcircuit name missing'/)
      go to 40
c
c  end
c
 5000 if (nsbckt.eq.0) go to 5010
      nsbckt=0
      write (iofile,5001)
 5001 format('0*error*:  .ends  card missing'/)
      nogo=1
 5010 call clrmem(ifield)
      call clrmem(icode)
      call clrmem(idelim)
      call clrmem(icolum)
      call clrmem(isbckt)
      if (nfour.eq.0) call clrmem(ifour)
      if (nsens.eq.0) call clrmem(isens)
 6000 call second(t2)
      rstats(1)=t2-t1
      return
      end
