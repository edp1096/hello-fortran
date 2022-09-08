      subroutine update(vinit,loct,node1,node2,nupda,icheck)
      implicit double precision (a-h,o-z)
c
c     this routine updates and limits the controlling variables for the
c nonlinear controlled sources.
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
      go to (40,10,40,20,30,50), initf
   10 vnew=vinit
      go to 70
   20 vnew=value(lx0+loct)
      go to 70
   30 vnew=value(lx1+loct)
      go to 70
   40 vnew=value(lvnim1+node1)-value(lvnim1+node2)
      go to 60
   50 call copy8(value(lx1+loct),value(lx0+loct),nupda)
      xfact=delta/delold(2)
      vnew=(1.0d0+xfact)*value(lx1+loct)-xfact*value(lx2+loct)
   60 if (dabs(vnew).le.1.0d0) go to 80
      delv=vnew-value(lx0+loct)
      if (dabs(delv).le.0.1d0) go to 80
      vlim=dmax1(dabs(0.1d0*value(lx0+loct)),0.1d0)
      vnew=value(lx0+loct)+dsign(dmin1(dabs(delv),vlim),delv)
      go to 70
   70 icheck=1
   80 value(lx0+loct)=vnew
      return
      end
