      subroutine dmpmat(anam)
      implicit double precision (a-h,o-z)
c
c      this routine dumps out the matrix and associated pointers.
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
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c spice version 2g.6  sccsid=cirdat 3/15/83
      common /cirdat/ locate(50),jelcnt(50),nunods,ncnods,numnod,nstop,
     1   nut,nlt,nxtrm,ndist,ntlin,ibr,numvs,numalt,numcyc
c spice version 2g.6  sccsid=mosarg 3/15/83
      common /mosarg/ vto,beta,gamma,phi,phib,cox,xnsub,xnfs,xd,xj,xld,
     1   xlamda,uo,uexp,vbp,utra,vmax,xneff,xl,xw,vbi,von,vdsat,qspof,
     2   beta0,beta1,cdrain,xqco,xqc,fnarrw,fshort,lev
c spice version 2g.6  sccsid=status 3/15/83
      common /status/ omega,time,delta,delold(7),ag(7),vt,xni,egfet,
     1   xmu,sfactr,mode,modedc,icalc,initf,method,iord,maxord,noncon,
     2   iterno,itemno,nosolv,modac,ipiv,ivmflg,ipostp,iscrch,iofile
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=knstnt 3/15/83
      common /knstnt/ twopi,xlog2,xlog10,root2,rad,boltz,charge,ctok,
     1   gmin,reltol,abstol,vntol,trtol,chgtol,eps0,epssil,epsox,
     2   pivtol,pivrel
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
c spice version 2g.6  sccsid=dc 3/15/83
      common /dc/ tcstar(2),tcstop(2),tcincr(2),icvflg,itcelm(2),kssop,
     1   kinel,kidin,kovar,kidout
c spice version 2g.6  sccsid=ac 3/15/83
      common /ac/ fstart,fstop,fincr,skw2,refprl,spw2,jacflg,idfreq,
     1   inoise,nosprt,nosout,nosin,idist,idprt
c spice version 2g.6  sccsid=tran 3/15/83
      common /tran/ tstep,tstop,tstart,delmax,tdmax,forfre,jtrflg
c spice version 2g.6  sccsid=outinf 3/15/83
      common /outinf/ xincr,string(15),xstart,yvar(8),itab(8),itype(8),
     1   ilogy(8),npoint,numout,kntr,numdgt
c spice version 2g.6  sccsid=cje 3/15/83
      common /cje/ maxtim,itime,icost
c spice version 2g.6  sccsid=debug 3/15/83
      common/debug/ idebug(20)
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      write (iofile,10) anam,mode,modedc,
     1 time,delta,icalc,iterno,initf,ipiv,iord,noncon,igoof,nogo
   10 format('0*debug*:  dmpmt called by ',a8,/,
     1   ' *debug*:  mode, mdc, time, delta, icalc, itr#, initf,',
     2             ' piv, ord, ncon, igoof, nogo =',/,
     3   ' *debug*:  ',2i5,1p2d10.2,8i5)
      call dmpmem(5hdmpmt)
c
c  dump out the *whole* thing
c
      call sizmem(irpt,irpts)
      write (iofile,16) nstop,nttbr,irpts
   16 format(' *debug*:  nstop, nttbr, size(irpt) = ',3i6,/,
     1   ' *debug*:   index  irpt  irow  jcol  jcpt       value',
     2          10x,' index  irpt  irow  jcol  jcpt       value')
      j=(irpts+1)/2
      istop=j
      do 30 i=1,istop
      j=j+1
      write (iofile,26)
     1   i,nodplc(irpt+i),nodplc(irowno+i),nodplc(jcolno+i),
     2   nodplc(jcpt+i),value(lvn+i),
     3   j,nodplc(irpt+j),nodplc(irowno+j),nodplc(jcolno+j),
     4   nodplc(jcpt+j),value(lvn+j)
   26 format(' *debug*:  ',5i6,1pd12.4,10x,5i6,1pd12.4)
   30 continue
cc 51 format(" *debug*:  irpt   = ",18i6)
cc    write (iofile,56) (nodplc(irowno+i),i=1,irpts)
cc 56 format(" *debug*:  irowno = ",18i6)
cc    write (iofile,61) (nodplc(jcolno+i),i=1,irpts)
cc 61 format(" *debug*:  jcolno = ",18i6)
cc    write (iofile,66) (nodplc(jcpt  +i),i=1,irpts)
cc 66 format(" *debug*:  jcpt   = ",18i6)
      write (iofile,71) (nodplc(irswpf+i),i=1,nstop)
   71 format(' *debug*:  irswpf = ',18i6)
      write (iofile,76) (nodplc(irswpr+i),i=1,nstop)
   76 format(' *debug*:  irswpr = ',18i6)
      write (iofile,81) (nodplc(icswpf+i),i=1,nstop)
   81 format(' *debug*:  icswpf = ',18i6)
      write (iofile,86) (nodplc(icswpr+i),i=1,nstop)
   86 format(' *debug*:  icswpr = ',18i6)
c
c
  500 return
      end
