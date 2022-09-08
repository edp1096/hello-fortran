      subroutine magphs(cvar,xmag,xphs)
      implicit double precision (a-h,o-z)
c
c     this routine computes the magnitude and phase of its complex arg-
c ument cvar, storing the results in xmag and xphs.
c
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
      complex cvar
c
c
      xreal=dble(real(cvar))
      ximag=dble(aimag(cvar))
      xmag=dsqrt(xreal*xreal+ximag*ximag)
      if (xmag.ge.1.0d-20) go to 10
      xmag=1.0d-20
      xphs=0.0d0
      return
   10 xphs=rad*datan2(ximag,xreal)
      return
      end
