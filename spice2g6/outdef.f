      subroutine outdef(ifld,mode,loct,ltype)
      implicit double precision (a-h,o-z)
c
c     this routine constructs the internal list element for an output
c variable defined on some input card.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      integer xxor
      dimension aout(19),aopts(5)
      data aout / 4hv   , 4hvm  , 4hvr  , 4hvi  , 4hvp  , 4hvdb ,
     1            4hi   , 4him  , 4hir  , 4hii  , 4hip  , 4hidb ,
     2            4honoi, 4hinoi, 4hhd2 , 4hhd3 , 4hdim2, 4hsim2,
     3            4hdim3 /
      data aopts / 1hm, 1hr, 1hi, 1hp, 1hd /
      data alprn, acomma, ablnk, aletv / 1h(, 1h,, 1h , 1hv /
c
      if (nodplc(icode+ifld).ne.1) go to 300
      anam=value(ifield+ifld)
      call move(anam,5,ablnk,1,4)
      do 10 i=1,19
      if (xxor(anam,aout(i)).ne.0) go to 10
      idout=i
      go to 20
   10 continue
      go to 300
c
c  further error checking
c
   20 if (mode.ge.3) go to 25
c...  dc or tran
      if ((idout.ne.1).and.(idout.ne.7)) go to 300
      go to 38
   25 if (mode.ge.4) go to 30
c...  ac
      if (idout.ge.13) go to 300
      go to 38
   30 if (mode.eq.5) go to 35
c...  noise
      if ((idout.ne.13).and.(idout.ne.14)) go to 300
      go to 38
c...  distortion
   35 if (idout.lt.15) go to 300
   38 ktype=0
      ltype=idout
      if (idout.lt.7) go to 40
      ktype=1
      ltype=ltype-6
      if (idout.lt.13) go to 40
      ktype=idout-11
      ltype=1
c
c  voltage output
c
   40 id=40+mode
      if (ktype.ne.0) go to 100
      if (nodplc(icode+ifld+1).ne.0) go to 300
      ifld=ifld+1
      n1=value(ifield+ifld)
      if (n1.lt.0) go to 300
      if(n1.gt.9999) go to 300
      n2=0
      adelim=value(idelim+ifld)
      if (adelim.eq.acomma) go to 45
      if (adelim.ne.ablnk) go to 50
   45 if (nodplc(icode+ifld+1).ne.0) go to 300
      ifld=ifld+1
      n2=value(ifield+ifld)
      if (n2.lt.0) go to 300
      if(n2.gt.9999) go to 300
   50 outnam=ablnk
      ipos=1
      call alfnum(n1,outnam,ipos)
      ipos=5
      call alfnum(n2,outnam,ipos)
      call find(outnam,id,loct,0)
      nodplc(loct+2)=n1
      nodplc(loct+3)=n2
      go to 400
c
c  current output
c
  100 if (ktype.ne.1) go to 200
      if (nodplc(icode+ifld+1).ne.1) go to 300
      ifld=ifld+1
      avsrc=value(ifield+ifld)
      achek=avsrc
      call move(achek,2,ablnk,1,7)
      if (achek.ne.aletv) go to 300
      call find(avsrc,id,loct,0)
      call find(avsrc,9,nodplc(loct+2),0)
      nodplc(loct+5)=1
      go to 400
c
c  noise or distortion outputs
c
  200 id=44
      if (ktype.ge.4) id=id+1
      if (value(idelim+ifld).ne.alprn) go to 220
      if (nodplc(icode+ifld+1).ne.1) go to 300
      ifld=ifld+1
      atype=value(ifield+ifld)
      call move(atype,2,ablnk,1,7)
      do 210 i=1,5
      if (atype.ne.aopts(i)) go to 210
      ltype=i+1
      go to 220
  210 continue
      go to 300
  220 call find(anam,id,loct,0)
      nodplc(loct+2)=0
      nodplc(loct+5)=ktype
      go to 400
c
c  errors
c
  300 igoof=1
c
c  finished
c
  400 return
      end
