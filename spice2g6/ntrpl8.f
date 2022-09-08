      subroutine ntrpl8(locx,locy,numpnt)
      implicit double precision (a-h,o-z)
c
c     this routine interpolates the analysis data to obtain the values
c printed and/or plotted, using linear interpolation.
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
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=outinf 3/15/83
      common /outinf/ xincr,string(15),xstart,yvar(8),itab(8),itype(8),
     1   ilogy(8),npoint,numout,kntr,numdgt
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
c  for dc transfer curve, no interpolation necessary
c
      if(mode.ne.1) go to 4
      numpnt=icalc
      loco=loutpt
      do 3 i=1,numpnt
      locyt=locy
      value(locx+i)=value(loco+1)
      do 2 k=1,kntr
      iseq=itab(k)
      iseq=nodplc(iseq+4)
      value(locyt+i)=value(loco+iseq)
      locyt=locyt+npoint
    2 continue
      loco=loco+numout
    3 continue
      return
    4 continue
      xvar=xstart
      xvtol=xincr*1.0d-5
      ippnt=0
      icpnt=2
      loco1=loutpt
      loco2=loco1+numout
      if (icalc.lt.2) go to 50
   10 x1=value(loco1+1)
      x2=value(loco2+1)
      dx1x2=x1-x2
   20 if (xincr.lt.0.0d0) go to 24
      if (xvar.le.(x2+xvtol)) go to 30
      go to 28
   24 if (xvar.ge.(x2+xvtol)) go to 30
   28 if (icpnt.ge.icalc) go to 100
      icpnt=icpnt+1
      loco1=loco2
      loco2=loco1+numout
      go to 10
   30 ippnt=ippnt+1
      value(locx+ippnt)=xvar
      dxx1=xvar-x1
      locyt=locy
      do 40 i=1,kntr
      iseq=itab(i)
      iseq=nodplc(iseq+4)
      v1=value(loco1+iseq)
      v2=value(loco2+iseq)
      yvr=v1+(v1-v2)*dxx1/dx1x2
      tol=dmin1(dabs(v1),dabs(v2))*1.0d-10
      if (dabs(yvr).le.tol) yvr=0.0d0
      value(locyt+ippnt)=yvr
      locyt=locyt+npoint
   40 continue
      if (ippnt.ge.npoint) go to 100
      xvar=xstart+dble(ippnt)*xincr
      if (dabs(xvar).ge.dabs(xvtol)) go to 20
      xvar=0.0d0
      go to 20
c
c  special handling if icalc = 1
c
c...  icalc=1;  just copy over the single point and return
   50 ippnt=1
      value(locx+ippnt)=xvar
      locyt=locy
      do 60 i=1,kntr
      iseq=itab(i)
      iseq=nodplc(iseq+4)
      value(locyt+ippnt)=value(loco1+iseq)
      locyt=locyt+npoint
   60 continue
      go to 100
c
c  return
c
  100 numpnt=ippnt
      return
      end
