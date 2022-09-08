      subroutine outnam(loc,ktype,string,ipos)
      implicit double precision (a-h,o-z)
c
c     this routine constructs the 'name' for the output variable indi-
c cated by loc, adding the characters to the character array 'string',
c beginning with the position marked by ipos.
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      dimension string(1)
      dimension aout(19),lenout(19),aopt(5),lenopt(5)
      data aout / 6hv     , 6hvm    , 6hvr    , 6hvi    , 6hvp    ,
     1            6hvdb   , 6hi     , 6him    , 6hir    , 6hii    ,
     2            6hip    , 6hidb   , 6honoise, 6hinoise, 6hhd2   ,
     1            6hhd3   , 6hdim2  , 6hsim2  , 6hdim3   /
      data lenout / 1,2,2,2,2,3,1,2,2,2,2,3,6,6,3,3,4,4,4 /
      data aopt / 5hmag  , 5hreal , 5himag , 5hphase, 5hdb    /
      data lenopt / 3,4,4,5,2 /
      data alprn, acomma, arprn, ablnk / 1h(, 1h,, 1h), 1h  /
c
c
      ioutyp=nodplc(loc+5)
      if (ioutyp.ge.2) go to 10
      lout=ktype+ioutyp*6
      go to 20
   10 lout=ioutyp+11
   20 call move(string,ipos,aout(lout),1,lenout(lout))
      ipos=ipos+lenout(lout)
      if (ioutyp.ge.2) go to 200
      call move(string,ipos,alprn,1,1)
      ipos=ipos+1
      if (ioutyp.ne.0) go to 100
      node1=nodplc(loc+2)
      call alfnum(nodplc(junode+node1),string,ipos)
      node2=nodplc(loc+3)
      if (node2.eq.1) go to 30
      call move(string,ipos,acomma,1,1)
      ipos=ipos+1
      call alfnum(nodplc(junode+node2),string,ipos)
   30 call move(string,ipos,arprn,1,1)
      ipos=ipos+1
      go to 1000
c
  100 locv=nodplc(loc+1)
      anam=value(locv)
      achar=ablnk
      do 110 i=1,8
      call move(achar,1,anam,i,1)
      if (achar.eq.ablnk) go to 120
      call move(string,ipos,achar,1,1)
      ipos=ipos+1
  110 continue
  120 call move(string,ipos,arprn,1,1)
      ipos=ipos+1
      go to 1000
c
  200 if (ktype.eq.1) go to 1000
      call move(string,ipos,alprn,1,1)
      ipos=ipos+1
      call move(string,ipos,aopt(ktype-1),1,lenopt(ktype-1))
      ipos=ipos+lenopt(ktype-1)
      call move(string,ipos,arprn,1,1)
      ipos=ipos+1
c
c  finished
c
 1000 return
      end
