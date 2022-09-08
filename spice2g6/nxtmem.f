      integer function nxtmem(memwds)
c
c.. function returns the in nxtmem the next available memory size
c.. (which must be evenly divisible by 'nwd4, nwd8, and nwd16' as
c.. defined in setmem
c
      nxtmem=((memwds+1999)/2000)*2000
      return
      end
