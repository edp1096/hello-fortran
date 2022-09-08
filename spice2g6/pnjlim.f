      subroutine pnjlim(vnew,vold,vt,vcrit,icheck)
      implicit double precision (a-h,o-z)
c
c     this routine limits the change-per-iteration of device pn-junction
c voltages.
c
      if (vnew.le.vcrit) go to 30
      vlim=vt+vt
      delv=vnew-vold
      if (dabs(delv).le.vlim) go to 30
      if (vold.le.0.0d0) go to 20
      arg=1.0d0+delv/vt
      if (arg.le.0.0d0) go to 10
      vnew=vold+vt*dlog(arg)
      go to 100
   10 vnew=vcrit
      go to 100
   20 vnew=vt*dlog(vnew/vt)
      go to 100
   30 icheck=0
c
c  finished
c
  100 return
      end
