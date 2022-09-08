      subroutine sorstp(itlim)
      implicit double precision (a-h,o-z)
c
c     this routine uses the source stepping method to solve the dc
c     operating point
c
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
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
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
      bound=1.0d0/64
      fractn=1.0d0/16
c
c  step down sources
c
   10 fractn=fractn*2.0d0
      sfactr=sfactr*fractn
      if (sfactr.lt.bound) go to 100
      initf=2
      call iter8(itlim)
      rstats(6)=rstats(6)+iterno
      if (igoof.ne.0) go to 10
      fractn=2.0d0
c
c  step up sources
c
   20 sfactr=sfactr*fractn
      if (sfactr.le.1.0d0) go to 30
      sfactr=1.0d0
   30 initf=3
      call iter8(itlim)
      rstats(6)=rstats(6)+iterno
      if ((igoof.eq.0).and.(sfactr.eq.1.0d0)) go to 200
      if (igoof.eq.0) go to 20
c
c  step down if step up failed
c
   40 fractn=dsqrt(fractn)
      if (fractn.lt.1.0001d0) go to 100
      sfactr=sfactr/fractn
      initf=3
      call iter8(itlim)
      rstats(6)=rstats(6)+iterno
      if (igoof.ne.0) go to 40
      go to 20
c
c   finish with source stepping method
c
  100 igoof=1
      write(iofile,110)
  110 format('0 source stepping method failed')
  200 initf=2
      return
      end
