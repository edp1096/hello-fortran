      subroutine shlsrt(a,n)
      implicit double precision (a-h,o-z)
c
c     this routine sorts the array a using a shell sort algorithm.
c
      dimension a(n)
      integer h
c
c
c...  compute best starting step size
      h=1
   10 h=3*h+1
      if (h.lt.n) go to 10
c...  back off two times
      h=(h-1)/3
      h=(h-1)/3
      h=max0(h,1)
c
c  shell sort
c
   20 j=h+1
      go to 60
   30 i=j-h
c...  ak = record key;  ar = record
      ak=a(j)
      ar=ak
   40 if (ak.ge.a(i)) go to 50
      a(i+h)=a(i)
      i=i-h
      if (i.ge.1) go to 40
   50 a(i+h)=ar
      j=j+1
   60 if (j.le.n) go to 30
      h=(h-1)/3
      if (h.ne.0) go to 20
      return
      end
