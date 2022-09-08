      subroutine trunc(delnew)
      implicit double precision (a-h,o-z)
c
c     this routine determines the new transient stepsize by either
c calling terr to estimate the local truncation error, or by checking
c on the number of iterations needed to converge at the last timepoint.
c
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
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      if (lvltim.ne.0) go to 5
      delnew=dmin1(tstep,delmax)
      return
    5 if (lvltim.ne.1) go to 10
      delnew=delta
      if (iterno.gt.itl3) return
      delnew=dmin1(2.0d0*delta,tstep,delmax)
      return
c
c  capacitors
c
   10 delnew=1.0d20
      loc=locate(2)
   20 if ((loc.eq.0).or.(nodplc(loc+12).ne.0)) go to 30
      loct=nodplc(loc+8)
      call terr(loct,delnew)
      loc=nodplc(loc)
      go to 20
c
c  inductors
c
   30 loc=locate(3)
   40 if ((loc.eq.0).or.(nodplc(loc+14).ne.0)) go to 50
      loct=nodplc(loc+11)
      call terr(loct,delnew)
      loc=nodplc(loc)
      go to 40
c
c  diodes
c
   50 loc=locate(11)
   60 if ((loc.eq.0).or.(nodplc(loc+16).ne.0)) go to 70
      loct=nodplc(loc+11)
      call terr(loct+3,delnew)
      loc=nodplc(loc)
      go to 60
c
c  bjts
c
   70 loc=locate(12)
   80 if ((loc.eq.0).or.(nodplc(loc+36).ne.0)) go to 90
      loct=nodplc(loc+22)
      call terr(loct+8,delnew)
      call terr(loct+10,delnew)
      call terr(loct+12,delnew)
      loc=nodplc(loc)
      go to 80
c
c  jfets
c
   90 loc=locate(13)
  100 if ((loc.eq.0).or.(nodplc(loc+25).ne.0)) go to 110
      loct=nodplc(loc+19)
      call terr(loct+9,delnew)
      call terr(loct+11,delnew)
      loc=nodplc(loc)
      go to 100
c
c  mosfets
c
  110 loc=locate(14)
  120 if ((loc.eq.0).or.(nodplc(loc+33).ne.0)) go to 200
      loct=nodplc(loc+26)
      call terr(loct+12,delnew)
      call terr(loct+14,delnew)
      call terr(loct+16,delnew)
      loc=nodplc(loc)
      go to 120
c
c  delta is allowed only to double at each timepoint
c
  200 delnew=dmin1(2.0d0*delta,delnew,delmax)
      return
      end
