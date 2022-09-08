      subroutine swapij(i1,i2,j1,j2)
      implicit double precision (a-h,o-z)
c spice version 2g.6  sccsid=tabinf 3/15/83
      common /tabinf/ ielmnt,isbckt,nsbckt,iunsat,nunsat,itemps,numtem,
     1   isens,nsens,ifour,nfour,ifield,icode,idelim,icolum,insize,
     2   junode,lsbkpt,numbkp,iorder,jmnode,iur,iuc,ilc,ilr,numoff,isr,
     3   nmoffc,iseq,iseq1,neqn,nodevs,ndiag,iswap,iequa,macins,lvnim1,
     4   lx0,lvn,lynl,lyu,lyl,lx1,lx2,lx3,lx4,lx5,lx6,lx7,ld0,ld1,ltd,
     5   imynl,imvn,lcvn,nsnod,nsmat,nsval,icnod,icmat,icval,
     6   loutpt,lpol,lzer,irswpf,irswpr,icswpf,icswpr,irpt,jcpt,
     7   irowno,jcolno,nttbr,nttar,lvntmp
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c     swap rows i1 and i2
c
      loc1=nodplc(jcpt+i1)
      loc2=nodplc(jcpt+i2)
      nodplc(jcpt+i1)=loc2
      nodplc(jcpt+i2)=loc1
c
c     check if end of row
c
    5 if (loc1.le.0.and.loc2.le.0) go to 80
c
c     check swap type
c
      if (loc1.eq.0) go to 20
      if (loc2.eq.0) go to 10
      if (nodplc(jcolno+loc1)-nodplc(jcolno+loc2)) 10,15,20
   10 ktype=-1
      j=nodplc(jcolno+loc1)
      go to 25
   15 ktype=0
      j=nodplc(jcolno+loc1)
      go to 25
   20 ktype=1
      j=nodplc(jcolno+loc2)
c
c     find pointer to entry (i1,j)
c
   25 loc=j
   30 lsav1=loc
      loc=nodplc(irpt+loc)
      if (loc.eq.0) go to 40
      if ((nodplc(irowno+loc)-i1).lt.0) go to 30
c
c     find pointer to entry (i2,j)
c
   40 loc=j
   45 lsav2=loc
      loc=nodplc(irpt+loc)
      if (loc.eq.0) go to 55
      if ((nodplc(irowno+loc)-i2).lt.0) go to 45
c
c     branch for col j in row i1, in both row i1 and i2, or in row i2
c
   55 if (ktype) 60,70,75
c
c     entry (i1,j)
c
   60 if (lsav1.eq.lsav2) go to 65
      loc=nodplc(irpt+lsav2)
      nodplc(irpt+lsav2)=loc1
      nodplc(irpt+lsav1)=nodplc(irpt+loc1)
      nodplc(irpt+loc1)=loc
   65 nodplc(irowno+loc1)=i2
      loc1=nodplc(jcpt+loc1)
      go to 5
c
c     entries (i1,j) and (i2,j)
c
   70 nodplc(irpt+lsav1)=loc2
      nodplc(irpt+lsav2)=loc1
      loc=nodplc(irpt+loc1)
      nodplc(irpt+loc1)=nodplc(irpt+loc2)
      nodplc(irpt+loc2)=loc
      nodplc(irowno+loc1)=i2
      nodplc(irowno+loc2)=i1
      loc1=nodplc(jcpt+loc1)
      loc2=nodplc(jcpt+loc2)
      go to 5
c
c     entry (i2,j)
c
   75 if (lsav1.eq.lsav2) go to 78
      loc=nodplc(irpt+lsav1)
      nodplc(irpt+lsav1)=loc2
      nodplc(irpt+lsav2)=nodplc(irpt+loc2)
      nodplc(irpt+loc2)=loc
   78 nodplc(irowno+loc2)=i1
      loc2=nodplc(jcpt+loc2)
      go to 5
c
c     swap columns j1 and j2
c
   80 loc1=nodplc(irpt+j1)
      loc2=nodplc(irpt+j2)
      nodplc(irpt+j1)=loc2
      nodplc(irpt+j2)=loc1
c
c     check for end of column
c
   85 if (loc1.le.0.and.loc2.le.0) go to 160
c
c     check swap type
c
      if (loc1.eq.0) go to 100
      if (loc2.eq.0) go to 90
      if (nodplc(irowno+loc1)-nodplc(irowno+loc2)) 90,95,100
   90 ktype=-1
      i=nodplc(irowno+loc1)
      go to 105
   95 ktype=0
      i=nodplc(irowno+loc1)
      go to 105
  100 ktype=1
      i=nodplc(irowno+loc2)
c
c     find pointer to entry (i,j1)
c
  105 loc=i
  110 lsav1=loc
      loc=nodplc(jcpt+loc)
      if (loc.eq.0) go to 120
      if ((nodplc(jcolno+loc)-j1).lt.0) go to 110
c
c     find pointer to entry (i,j2)
c
  120 loc=i
  125 lsav2=loc
      loc=nodplc(jcpt+loc)
      if(loc.eq.0) go to 135
      if ((nodplc(jcolno+loc)-j2).lt.0) go to 125
c
c     branch for row i in col j1, in both col"s j1 and j2, or in col j2
c
  135 if (ktype) 140,150,155
c
c     entry (i,j1)
c
  140 if (lsav1.eq.lsav2) go to 145
      loc=nodplc(jcpt+lsav2)
      nodplc(jcpt+lsav2)=loc1
      nodplc(jcpt+lsav1)=nodplc(jcpt+loc1)
      nodplc(jcpt+loc1)=loc
  145 nodplc(jcolno+loc1)=j2
      loc1=nodplc(irpt+loc1)
      go to 85
c
c     entries (i1,j) and (i2,j)
c
  150 nodplc(jcpt+lsav1)=loc2
      nodplc(jcpt+lsav2)=loc1
      loc=nodplc(jcpt+loc1)
      nodplc(jcpt+loc1)=nodplc(jcpt+loc2)
      nodplc(jcpt+loc2)=loc
      nodplc(jcolno+loc1)=j2
      nodplc(jcolno+loc2)=j1
      loc1=nodplc(irpt+loc1)
      loc2=nodplc(irpt+loc2)
      go to 85
c
c     entry (i,j2)
c
  155 if (lsav1.eq.lsav2) go to 158
      loc=nodplc(jcpt+lsav1)
      nodplc(jcpt+lsav1)=loc2
      nodplc(jcpt+lsav2)=nodplc(jcpt+loc2)
      nodplc(jcpt+loc2)=loc
  158 nodplc(jcolno+loc2)=j1
      loc2=nodplc(irpt+loc2)
      go to 85
  160 return
      end
