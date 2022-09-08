      subroutine nxtpwr(pwrseq,pdim)
      implicit double precision (a-h,o-z)
c
c     this routine determines the 'next' set of exponents for the
c different dimensions of a polynomial.
c
      integer pwrseq(1),pdim,psum
c
c
      if (pdim.eq.1) go to 80
      k=pdim
   10 if (pwrseq(k).ne.0) go to 20
      k=k-1
      if (k.ne.0) go to 10
      go to 80
   20 if (k.eq.pdim) go to 30
      pwrseq(k)=pwrseq(k)-1
      pwrseq(k+1)=pwrseq(k+1)+1
      go to 100
   30 km1=k-1
      do 40 i=1,km1
      if (pwrseq(i).ne.0) go to 50
   40 continue
      pwrseq(1)=pwrseq(pdim)+1
      pwrseq(pdim)=0
      go to 100
   50 psum=1
      k=pdim
   60 if (pwrseq(k-1).ge.1) go to 70
      psum=psum+pwrseq(k)
      pwrseq(k)=0
      k=k-1
      go to 60
   70 pwrseq(k)=pwrseq(k)+psum
      pwrseq(k-1)=pwrseq(k-1)-1
      go to 100
   80 pwrseq(1)=pwrseq(1)+1
c
c  finished
c
  100 return
      end
