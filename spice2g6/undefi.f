      subroutine undefi(array,length)
      implicit double precision (a-h,o-z)
c
      dimension array(1)
c     this routine undefines the memory locations indicated by array(1)
c through array(length).
c
      data aundef /2h.u/
      if (length.eq.0) return
      do 10 i=1,length
      array(i)=aundef
   10 continue
      return
      end
