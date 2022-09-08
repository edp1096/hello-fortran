      subroutine getcje
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
      call second(xtime)
      itime=xtime
      return
      end
