      subroutine cmult(xr,xi,yr,yi,cr,ci)
c.. ok if cr and ci are really xr and xi or yr and yi
      implicit double precision (a-h,o-z)
      xrtemp=xr
      xitemp=xi
      yrtemp=yr
      yitemp=yi
      cr=xrtemp*yrtemp-xitemp*yitemp
      ci=xitemp*yrtemp+xrtemp*yitemp
      return
      end
