      subroutine second(t1)
      implicit double precision (a-h,o-z)
      dimension ibuff(4)
      call times (ibuff)
      t1=dble (ibuff(1))/60.0d0
      return
      end
