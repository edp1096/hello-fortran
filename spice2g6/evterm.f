      subroutine evterm(val,arg,iexp)
      implicit double precision (a-h,o-z)
c
c     this routine evaluates one term of a polynomial.
c
      jexp=iexp+1
      if (jexp.ge.6) go to 60
      go to (10,20,30,40,50), jexp
   10 val=1.0d0
      go to 100
   20 val=arg
      go to 100
   30 val=arg*arg
      go to 100
   40 val=arg*arg*arg
      go to 100
   50 val=arg*arg
      val=val*val
      go to 100
   60 if (arg.eq.0.0d0) go to 70
      argexp=dble(iexp)*dlog(dabs(arg))
      if (argexp.lt.-200.0d0) go to 70
      val=dexp(argexp)
      if((iexp/2)*2.eq.iexp) go to 100
      val=dsign(val,arg)
      go to 100
   70 val=0.0d0
c
c  finished
c
  100 return
      end
