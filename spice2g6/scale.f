      subroutine scale(xmin,xmax,n,xminp,xmaxp,del)
      implicit double precision (a-h,o-z)
c
c     this routine determines the 'optimal' scale to use for the plot of
c some output variable.
c
c
c  adapted from algorithm 463 of 'collected algorithms of the cacm'
c
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
      integer xxor
      dimension vint(5)
      data vint / 1.0d0,2.0d0,5.0d0,10.0d0,20.0d0 /
      data eps / 1.0d-12 /
c
c
c...  trap too-small data spread
c***********************************************************
c  temporily check 'equality' this way
      if(xmin.eq.0.0d0.and.xmax.eq.0.0d0) go to 4
      if(dabs((xmax-xmin)/dmax1(dabs(xmin),dabs(xmax))).ge.1.0d-4)
     1  go to 10
    4 continue
      if (xmin.ge.0.0d0) go to 5
      xmax=0.5d0*xmin+eps
      xmin=1.5d0*xmin-eps
      go to 10
    5 xmax=1.5d0*xmin+eps
      xmin=0.5d0*xmin-eps
c...  find approximate interval size, normalized to [1,10]
   10 a=(xmax-xmin)/dble(n)
      nal=idint(dlog10(a))
      if (a.lt.1.0d0) nal=nal-1
      xfact=dexp(xlog10*dble(nal))
      b=a/xfact
c...  find closest permissible interval size
      do 20 i=1,3
      if (b.lt.(vint(i)+eps)) go to 30
   20 continue
      i=4
c...  compute interval size
   30 del=vint(i)*xfact
      fm1=xmin/del
      m1=fm1
      if (fm1.lt.0.0d0) m1=m1-1
      if (dabs(dble(m1)+1.0d0-fm1).lt.eps) m1=m1+1
c...  compute new maximum and minimum limits
      xminp=del*dble(m1)
      fm2=xmax/del
      m2=fm2+1.0d0
      if (fm2.lt.(-1.0d0)) m2=m2-1
      if (dabs(fm2+1.0d0-dble(m2)).lt.eps) m2=m2-1
      xmaxp=del*dble(m2)
      np=m2-m1
c...  check whether another loop required
      if (np.le.n) go to 40
      i=i+1
      go to 30
c...  do final adjustments and correct for roundoff error(s)
   40 nx=(n-np)/2
      xminp=dmin1(xmin,xminp-dble(nx)*del)
      xmaxp=dmax1(xmax,xminp+dble(n)*del)
      return
      end
