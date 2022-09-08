      integer function locf(ivar)
      implicit double precision (a-h,o-z)
      external loc
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
      dimension ivar(1)
      iabsa=loc(ivar)
      locf=iabsa/4
      if (iabsa.eq.4*locf) return
      write(iofile,100) iabsa
  100 format ('0*error*: system error, address ',i10,
     1   ' is not on 4-byte boundary')
      return
      end
