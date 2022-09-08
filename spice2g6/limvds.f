      subroutine limvds(vnew,vold)
      implicit double precision (a-h,o-z)
c
c     this routine limits the per-iteration change of fet vds.
c
      if (vold.lt.3.5d0) go to 200
c
      if (vnew.le.vold) go to 100
      vnew=dmin1(vnew,3.0d0*vold+2.0d0)
      go to 500
  100 if (vnew.lt.3.5d0) vnew=dmax1(vnew,2.0d0)
      go to 500
c
  200 if (vnew.le.vold) go to 300
      vnew=dmin1(vnew,4.0d0)
      go to 500
  300 vnew=dmax1(vnew,-0.5d0)
c
  500 return
      end
