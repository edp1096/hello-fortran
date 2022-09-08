      subroutine getlin
      implicit double precision (a-h,o-z)
c
c     this routine reads the next line of input into the array afield.
c if end-of-file is found, the variable keof is set to 1.
c
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c
c
      call copy8(afield,oldlin,15)
      read(5,6,end=10) (afield(i),i=1,10)
      go to 100
    6 format(10a8)
   10 keof=1
  100 call ushift(afield)
      return
      end
