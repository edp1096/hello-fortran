      subroutine fetlim(vnew,vold,vto)
      implicit double precision (a-h,o-z)
c
c     this routine limits the per-iteration change of fet voltages.
c
c
c      three regions of operation are identified:
c
c                  v < vto        definitely off
c            vto <= v <= vto+3.5d0    off or on depending on vbs
c        vto+3.5d0 < v              definitely on
c
      vtsthi=dabs(2.0d0*(vold-vto))+2.0d0
      vtstlo=vtsthi/2.0d0+2.0d0
      vtox=vto+3.5d0
      delv=vnew-vold
c
      if (vold.lt.vto) go to 300
      if (vold.lt.vtox) go to 200
c
c  on ...
c
      if (delv.gt.0.0d0) go to 120
c...  going off
      if (vnew.lt.vtox) go to 110
      if (-delv.le.vtstlo) go to 500
      vnew=vold-vtstlo
      go to 500
  110 vnew=dmax1(vnew,vto+2.0d0)
      go to 500
c...  staying on
  120 if (delv.lt.vtsthi) go to 500
      vnew=vold+vtsthi
      go to 500
c
c  middle region ...
c
  200 if (delv.gt.0.0d0) go to 210
c...  decreasing
      vnew=dmax1(vnew,vto-0.5d0)
      go to 500
c...  increasing
  210 vnew=dmin1(vnew,vto+4.0d0)
      go to 500
c
c  off ...
c
  300 if (delv.gt.0.0d0) go to 310
      if (-delv.le.vtsthi) go to 500
      vnew=vold-vtsthi
      go to 500
  310 vtemp=vto+0.5d0
      if (vnew.gt.vtemp) go to 320
      if (delv.le.vtstlo) go to 500
      vnew=vold+vtstlo
      go to 500
  320 vnew=vtemp
c
c  finished
c
  500 return
      end
