      subroutine plot(numpnt,locx,locy,locv)
      implicit double precision (a-h,o-z)
c
c     this routine generates the line-printer plots.
c
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
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
      integer xxor
      dimension ycoor(5,8),icoor(8),delplt(8)
      dimension agraph(13),aplot(13)
      dimension asym(2),pmin(8),jcoor(8)
      data ablnk, aletx, aper / 1h , 1hx, 1h. /
      data asym1, asym2, arprn / 8h(-------, 8h--------, 1h) /
      data pltsym / 8h*+=$0<>? /
c
c
      iwide=1
      nwide=101
      nwide4=25
      if(lwidth.gt.80) go to 3
      iwide=0
      nwide=57
      nwide4=14
    3 if (numpnt.le.0) go to 400
      do 5 i=1,13
      agraph(i)=ablnk
    5 continue
      do 7 i=1,5
      ispot=1+nwide4*(i-1)
      call move(agraph,ispot,aper,1,1)
    7 continue
      locyt=locy
      lspot=locv-1
      mltscl=0
      if (value(locv).eq.0.0d0) mltscl=1
      do 235 k=1,kntr
      lspot=lspot+2
      ymin=value(lspot)
      ymax=value(lspot+1)
      if (ymin.ne.0.0d0) go to 10
      if (ymax.ne.0.0d0) go to 10
      go to 100
   10 ymin1=dmin1(ymin,ymax)
      ymax1=dmax1(ymin,ymax)
   30 if (ilogy(k).eq.1) go to 40
      ymin1=dlog10(dmax1(ymin1,1.0d-20))
      ymax1=dlog10(dmax1(ymax1,1.0d-20))
      del=dmax1(ymax1-ymin1,0.0001d0)/4.0d0
      go to 50
   40 del=dmax1(ymax1-ymin1,1.0d-20)/4.0d0
   50 ymin=ymin1
      ymax=ymax1
      go to 200
c
c  determine max and min values
c
  100 ymax1=value(locyt+1)
      ymin1=ymax1
      if (numpnt.eq.1) go to 150
      do 110 i=2,numpnt
      ymin1=dmin1(ymin1,value(locyt+i))
      ymax1=dmax1(ymax1,value(locyt+i))
  110 continue
c
c  scaling
c
  150 call scale(ymin1,ymax1,4,ymin,ymax,del)
c
c  determine coordinates
c
  200 ycoor(1,k)=ymin
      pmin(k)=ymin
      small=del*1.0d-4
      if (dabs(ycoor(1,k)).le.small) ycoor(1,k)=0.0d0
      do 210 i=1,4
      ycoor(i+1,k)=ycoor(i,k)+del
      if (dabs(ycoor(i+1,k)).le.small) ycoor(i+1,k)=0.0d0
  210 continue
      if (ilogy(k).eq.1) go to 230
      do 220 i=1,5
  220 ycoor(i,k)=dexp(xlog10*ycoor(i,k))
  230 delplt(k)=del/dble(nwide4)
      locyt=locyt+npoint
  235 continue
c
c  count distinct coordinates
c
      icoor(1)=1
      jcoor(1)=1
      numcor=1
      if (kntr.eq.1) go to 290
      do 250 i=2,kntr
      do 245 j=1,numcor
      l=jcoor(j)
c...  coordinates are *equal* if the most significant 24 bits agree
      do 240 k=1,5
c*****************************************************************
c  temporarily check 'equality' this way
      y1=ycoor(k,i)
      y2=ycoor(k,l)
      if(y1.eq.0.0d0.and.y2.eq.0.0d0) go to 240
      if(dabs((y1-y2)/dmax1(dabs(y1),dabs(y2))).ge.1.0d-7) go to 245
  240 continue
      icoor(i)=l
      go to 250
  245 continue
      icoor(i)=i
      numcor=numcor+1
      jcoor(numcor)=i
  250 continue
c
c  print coordinates
c
  260 do 280 i=1,numcor
      asym(1)=asym1
      asym(2)=asym2
      ipos=2
      do 270 j=1,kntr
      if (icoor(j).ne.jcoor(i)) go to 270
      call move(asym,ipos,pltsym,j,1)
      ipos=ipos+1
  270 continue
      call move(asym,ipos,arprn,1,1)
      k=jcoor(i)
      if(iwide.ne.0) write(iofile,271) asym,(ycoor(j,k),j=1,5)
  271 format(/2a8,4h----,1pd12.3,4(15x,d10.3)/26x,51(2h -))
      if(iwide.eq.0) write(iofile,273) asym,(ycoor(j,k),j=1,5)
  273 format(/2a8,1x,1pd10.3,3(4x,d10.3),1x,d10.3/22x,29(2h -))
  280 continue
      go to 300
  290 if(iwide.ne.0) write(iofile,291) (ycoor(j,1),j=1,5)
  291 format(/20x,1pd12.3,4(15x,d10.3)/26x,51(2h -))
      if(iwide.eq.0) write(iofile,293) (ycoor(j,1),j=1,5)
  293 format(/15x,1pd12.3,3(4x,d10.3),1x,d10.3/22x,29(2h -))
c
c  plotting
c
  300 aspot=ablnk
      do 320 i=1,numpnt
      xvar=value(locx+i)
      locyt=locy
      call copy8(agraph,aplot,13)
      do 310 k=1,kntr
      yvr=value(locyt+i)
      ktmp=icoor(k)
      ymin1=pmin(ktmp)
      jpoint=idint((yvr-ymin1)/delplt(k)+0.5d0)+1
      if (jpoint.le.0) go to 306
      if (jpoint.gt.nwide) go to 306
      call move(aspot,1,aplot,jpoint,1)
      if (aspot.eq.ablnk) go to 303
      if (aspot.eq.aper) go to 303
      call move(aplot,jpoint,aletx,1,1)
      go to 306
  303 call move(aplot,jpoint,pltsym,k,1)
  306 locyt=locyt+npoint
  310 continue
      yvr=value(locy+i)
      if (ilogy(1).eq.1) go to 315
      yvr=dexp(xlog10*yvr)
  315 if(iwide.ne.0) write(iofile,316) xvar,yvr,aplot
  316 format(1pd10.3,3x,d10.3,3x,13a8)
      if(iwide.eq.0) write(iofile,317) xvar,yvr,(aplot(k),k=1,8)
  317 format(1pd10.3,1x,d10.3,1x,7a8,a1)
  320 continue
c
c  finished
c
      if(iwide.ne.0) write(iofile,331)
  331 format(26x,51(2h -)//)
      if(iwide.eq.0) write(iofile,332)
  332 format(22x,29(2h -)//)
      go to 500
c
c  too few points
c
  400 write (iofile,401)
  401 format('0warning:  too few points for plotting'/)
  500 write (iofile,501)
  501 format(1hy)
      return
      end
