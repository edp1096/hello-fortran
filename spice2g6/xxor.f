      integer function xxor(a,b)
      implicit double precision (a-h,o-z)
c
c     this routine computes a single-precision integer result which is
c the result of exclusive-or*ing the two real-valued arguments a and b
c together.
c
      xxor=1
      if(a.eq.b) xxor=0
      return
      end
