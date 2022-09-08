      integer function nxtevn(n)
c
c.. function returns the smallest value nxtevn greater than or equal to
c.. n which is evenly divisible by 'nwd4, nwd8, and nwd16' as defined
c.. in setmem
c
      nxtevn=((n+3)/4)*4
      return
      end
