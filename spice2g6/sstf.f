      subroutine sstf
      implicit double precision (a-h,o-z)
c
c     this routine computes the value of the small-signal transfer
c function specified by the user.
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
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      dimension string(5),save(3)
      data aslash, ablnk / 1h/, 1h  /
c
c  setup current vector for input resistance and transfer function
c
      call zero8(value(lvn+1),nstop)
      if (kidin.eq.10) go to 5
c...  voltage source input
      iptri=nodplc(kinel+6)
      value(lvn+iptri)=+1.0d0
      go to 20
c...  current source input
    5 noposi=nodplc(kinel+2)
      nonegi=nodplc(kinel+3)
      value(lvn+noposi)=-1.0d0
      value(lvn+nonegi)=+1.0d0
c
c  lu decompose and solve the system of circuit equations
c
c...  reorder the right-hand side
   20 call dcdcmp
      call dcsol
      value(lvn+1)=0.0d0
      do 25 i=1,nstop
      j=nodplc(icswpr+i)
      k=nodplc(irswpf+j)
      value(lvntmp+i)=value(lvn+k)
   25 continue
      call copy8(value(lvntmp+1),value(lvn+1),nstop)
c
c  evaluate transfer function
c
      if (nodplc(kovar+5).ne.0) go to 30
c...  voltage output
      noposo=nodplc(kovar+2)
      nonego=nodplc(kovar+3)
      trfn=value(lvn+noposo)-value(lvn+nonego)
      go to 40
c...  current output (through voltage source)
   30 iptro=nodplc(kovar+2)
      iptro=nodplc(iptro+6)
      trfn=value(lvn+iptro)
c
c  evaluate input resistance
c
   40 if (kidin.eq.9) go to 50
c...  current source input
      zin=value(lvn+nonegi)-value(lvn+noposi)
      go to 70
c...  voltage source input
   50 creal=value(lvn+iptri)
      if (dabs(creal).ge.1.0d-20) go to 60
      zin=1.0d20
      go to 70
   60 zin=-1.0d0/creal
c
c  setup current vector for output resistance
c
   70 call zero8(value(lvn+1),nstop)
      if (nodplc(kovar+5).ne.0) go to 80
c...  voltage output
      value(lvn+noposo)=-1.0d0
      value(lvn+nonego)=+1.0d0
      go to 90
   80 if (nodplc(kovar+2).ne.kinel) go to 85
      zout=zin
      go to 200
c...  current output (through voltage source)
   85 value(lvn+iptro)=+1.0d0
c
c  perform new forward and backward substitution
c
c...  reorder the right-hand side
   90 call dcsol
      value(lvn+1)=0.0d0
      do 95 i=1,nstop
      j=nodplc(icswpr+i)
      k=nodplc(irswpf+j)
      value(lvntmp+i)=value(lvn+k)
   95 continue
      call copy8(value(lvntmp+1),value(lvn+1),nstop)
c
c  evaluate output resistance
c
  100 if (nodplc(kovar+5).ne.0) go to 110
c...  voltage output
      zout=value(lvn+nonego)-value(lvn+noposo)
      go to 200
c...  current output (through voltage source)
  110 creal=value(lvn+iptro)
      if (dabs(creal).ge.1.0d-20) go to 120
      zout=1.0d20
      go to 200
  120 zout=-1.0d0/creal
c
c  print results
c
  200 do 210 i=1,5
      string(i)=ablnk
  210 continue
      ipos=1
      call outnam(kovar,1,string,ipos)
      call copy8(string,save,3)
      call move(string,ipos,aslash,1,1)
      ipos=ipos+1
      locv=nodplc(kinel+1)
      anam=value(locv)
      call move(string,ipos,anam,1,8)
      write (iofile,231) string,trfn,anam,zin,save,zout
  231 format(////,'0****     small-signal characteristics'//,
     1   1h0,5x,5a8,3h = ,1pd10.3,/,
     2   1h0,5x,'input resistance at ',a8,12x,3h = ,d10.3,/,
     3   1h0,5x,'output resistance at ',2a8,a3,3h = ,d10.3)
      return
      end
