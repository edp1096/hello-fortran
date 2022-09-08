      subroutine pheadr(aheadr)
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
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=miscel 3/15/83
      common /miscel/ atime,aprog(3),adate,atitle(10),defl,defw,defad,
     1  defas,rstats(50),iwidth,lwidth,nopage
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
c int3 (not used) is strictly for alignment.  f77 on unix craps out.
      integer int2,int3,nodpl2(128)
      equivalence (value(1),nodpl2(1))
      equivalence (value(1),nodplc(1),cvalue(1))
      dimension aheadr(10)
c
c  put out the header records onto the post-processing file
c  routine is used for all analysis modes (mode=1,2,3)
c
      dimension xtype(2)
      data xtype /4htime,4hfreq/
      data ablnk,aletv,aleti /1h ,1hv,1hi/
c
c file structure for post-processor
c
c record 1  title card (80 bytes), date (8 bytes), time (8 bytes) total-96 bytes
c record 2  number of output variables (including "sweep" variable)
c record 3  integer '4' (2 bytes)
c record 4  names of each output variable (8 bytes ea.)
c record 5  type of each output       0-no type
c                                     1-time
c                                     2-frequency
c                                     3-voltage
c                                     4-current
c                                     5-output noise
c                                     6-input noise
c                                     7-hd2    |
c                                     8-hd3    |
c                                     9-dim2   }   distortion outputs
c                                    10-sim2   |
c                                    11-dim3   |
c record 6  the location of each variable within each sweep point.
c           (normally just 1,2,3,4,... but needed if outputs are mixed up)
c record 6a 24 characters that are the plot title if record 3 is a '4'.
c record 7  output at first sweep point
c record 8  output at second sweep point
c record 9  .
c           .
c           .
c last record
c
c
      call getm8(ibuff,12)
      call copy8(aheadr(1),value(ibuff+1),10)
      value(ibuff+11)=adate
      value(ibuff+12)=atime
      call fwrite(value(ibuff+1),48)
      numout=nunods+jelcnt(9)
c force nused to be allocated by useless usage.
      int2 = numout
      int3 = numout
      info=4
      call getm8(inames,numout)
      call getm4(itypes,numout)
      call getm4(iseqs,numout)
      itype2=itypes*2
      iseq2=iseqs*2
      iknt=1
      nodpl2(iseq2+1)=1
c
c dc transfer curve (mode = 1):
c
      if(mode.ne.1) go to 10
      loc=itcelm(1)
      locv=nodplc(loc+1)
      value(inames+1)=value(locv)
      anam=ablnk
      call move(anam,1,value(locv),1,1)
      ityp=0
c voltage transfer becomes type 3 and current transfer becomes 4.
      if(anam.eq.aletv) ityp=3
      if(anam.eq.aleti) ityp=4
      nodpl2(itype2+1)=ityp
      go to 20
   10 value(inames+1)=xtype(mode-1)
      nodpl2(itype2+1)=mode-1
   20 do 30 i=2,nunods
      nodpl2(itype2+i)=3
      nodpl2(iseq2+i)=i
      value(inames+i)=ablnk
      ipos=1
      call alfnum(nodplc(junode+i),value(inames+i),ipos)
   30 continue
      loc=locate(9)
      iknt=nunods
   40 if(loc.eq.0) go to 50
      iknt=iknt+1
      nodpl2(itype2+iknt)=4
      nodpl2(iseq2+iknt)=iknt
      locv=nodplc(loc+1)
      value(inames+iknt)=value(locv)
      loc=nodplc(loc)
      go to 40
   50 int2=numout
      call fwrite(int2,1)
      int2=info
      call fwrite(int2,1)
      nwds=numout*4
      call fwrite(value(inames+1),nwds)
      call fwrite(nodpl2(itype2+1),numout)
      call fwrite(nodpl2(iseq2+1),numout)
      call fwrite(aprog(1),12)
      call clrmem(ibuff)
      call clrmem(inames)
      call clrmem(itypes)
      call clrmem(iseqs)
      return
      end
