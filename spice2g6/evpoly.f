      subroutine evpoly(result,itype,lcoef,ncoef,larg,
     1  narg,lexp)
      implicit double precision (a-h,o-z)
c
c     this routine evaluates a polynomial.  lcoef points to the coef-
c ficients, and larg points to the values of the polynomial argument(s).
c
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      if (itype) 100,200,300
c
c  integration (polynomial *must* be one-dimensional)
c
  100 result=0.0d0
      arg=1.0d0
      arg1=value(larg+1)
      do 110 i=1,ncoef
      arg=arg*arg1
      result=result+value(lcoef+i)*arg/dble(i)
  110 continue
      go to 1000
c
c  evaluation of the polynomial
c
  200 result=value(lcoef+1)
      if (ncoef.eq.1) go to 1000
      call zero4(nodplc(lexp+1),narg)
      do 220 i=2,ncoef
      call nxtpwr(nodplc(lexp+1),narg)
      if (value(lcoef+i).eq.0.0d0) go to 220
      arg=1.0d0
      do 210 j=1,narg
      call evterm(val,value(larg+j),nodplc(lexp+j))
      arg=arg*val
  210 continue
      result=result+value(lcoef+i)*arg
  220 continue
      go to 1000
c
c  partial derivative with respect to the itype*th variable
c
  300 result=0.0d0
      if (ncoef.eq.1) go to 1000
      call zero4(nodplc(lexp+1),narg)
      do 330 i=2,ncoef
      call nxtpwr(nodplc(lexp+1),narg)
      if (nodplc(lexp+itype).eq.0) go to 330
      if (value(lcoef+i).eq.0.0d0) go to 330
      arg=1.0d0
      do 320 j=1,narg
      if (j.eq.itype) go to 310
      call evterm(val,value(larg+j),nodplc(lexp+j))
      arg=arg*val
      go to 320
  310 call evterm(val,value(larg+j),nodplc(lexp+j)-1)
      arg=arg*dble(nodplc(lexp+j))*val
  320 continue
      result=result+value(lcoef+i)*arg
  330 continue
c
c  finished
c
 1000 return
      end
