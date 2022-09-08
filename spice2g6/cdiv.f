      subroutine cdiv(xr,xi,yr,yi,cr,ci)
c.. ok if cr and ci are really xr and xi or yr and yi
      implicit double precision (a-h,o-z)
      xrtemp=xr
      xitemp=xi
      yrtemp=yr
      yitemp=yi
      amag2=yrtemp*yrtemp+yitemp*yitemp
      cr=(xrtemp*yrtemp+xitemp*yitemp)/amag2
      ci=(xitemp*yrtemp-xrtemp*yitemp)/amag2
      return
      end
