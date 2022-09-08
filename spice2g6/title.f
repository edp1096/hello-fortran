      subroutine title(ifold,len,icom,coment)
      implicit double precision (a-h,o-z)
c
c     this routine writes a title on the output file.  ifold indicates
c whether the page eject should be to the next concave, convex, or any
c page fold depending on whether its value is <0, >0, or =0.  the page
c eject is suppressed (as is much of the heading) if the variable nopage
c is nonzero.
c
c spice version 2g.6  sccsid=tabinf 3/15/83
      common /tabinf/ ielmnt,isbckt,nsbckt,iunsat,nunsat,itemps,numtem,
     1   isens,nsens,ifour,nfour,ifield,icode,idelim,icolum,insize,
     2   junode,lsbkpt,numbkp,iorder,jmnode,iur,iuc,ilc,ilr,numoff,isr,
     3   nmoffc,iseq,iseq1,neqn,nodevs,ndiag,iswap,iequa,macins,lvnim1,
     4   lx0,lvn,lynl,lyu,lyl,lx1,lx2,lx3,lx4,lx5,lx6,lx7,ld0,ld1,ltd,
     5   imynl,imvn,lcvn,nsnod,nsmat,nsval,icnod,icmat,icval,
     6   loutpt,lpol,lzer,irswpf,irswpr,icswpf,icswpr,irpt,jcpt,
     7   irowno,jcolno,nttbr,nttar,lvntmp
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension coment(4)
c
c
      if(nopage.eq.1) go to 150
c
   30 if (len.le.80) go to 100
      write (iofile,31) adate,aprog,atime,(atitle(i),i=1,10)
   31 format(1h1,16(1h*),a8,1x,24(1h*),3a8,24(1h*),a8,16(1h*),//1h0,
     1   15a8/)
      if (icom.eq.0) go to 40
      write (iofile,36) coment,value(itemps+itemno)
   36 format(5h0****,17x,4a8,21x,'temperature =',f9.3,' deg c'/)
   40 write (iofile,41)
   41 format(1h0,121(1h*)//)
      go to 200
c
c
  100 write (iofile,101) adate,aprog,atime,(atitle(i),i=1,10)
  101 format(1h1,7(1h*),a8,1x,8(1h*),3a8,8(1h*),a8,5(1h*)//1h0,10a8/)
      if (icom.eq.0) go to 110
      write (iofile,106) coment,value(itemps+itemno)
  106 format(10h0****     ,4a8,' temperature =',f9.3,' deg c'/)
  110 write (iofile,111)
  111 format(1h0,71(1h*)//)
      go to 200
c
c
  150 if (icom.eq.0) go to 160
      write (iofile,106) coment,value(itemps+itemno)
      go to 200
  160 write (iofile,161) aprog
  161 format(1h0,3a8,/)
c
c  finished
c
  200 return
      end
