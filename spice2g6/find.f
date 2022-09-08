      subroutine find(aname,id,loc,iforce)
      implicit double precision (a-h,o-z)
c
c     this routine searches the list with number 'id' for an element
c with name 'aname'.  loc is set to point to the element.  if iforce is
c nonzero, then find expects to have to add the element to the list, and
c reports a fatal error if the element is found.  if subcircuit defini-
c tion is in progress (nonzero value for nsbckt), then find searches the
c current subcircuit definition list rather than the nominal element
c list.
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
c spice version 2g.6  sccsid=flags 3/15/83
      common /flags/ iprnta,iprntl,iprntm,iprntn,iprnto,limtim,limpts,
     1   lvlcod,lvltim,itl1,itl2,itl3,itl4,itl5,itl6,igoof,nogo,keof
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
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
c  index to the contents of the various lists:
c
c        list      contents
c        ----      --------
c
c          1       resistors
c          2       nonlinear capacitors
c          3       nonlinear inductors
c          4       mutual inductors
c          5       nonlinear voltage controlled current sources
c          6       nonlinear voltage controlled voltage sources
c          7       nonlinear current controlled current sources
c          8       nonlinear current controlled voltage sources
c          9       independent voltage sources
c         10       independent current sources
c         11       diodes
c         12       bipolar junction transistors
c         13       junction field-effect transistors (jfets)
c         14       metal-oxide-semiconductor junction fets (mosfets)
c         15       s-parameter 2-port network
c         16       y-parameter 2-port network
c         17       transmission lines
c         18       used for temperature sweeping
c         19       subcircuit calls
c         20       subcircuit definitions
c         21       diode model
c         22       bjt model
c         23       jfet model
c         24       mosfet model
c      25-30       <unused>
c         31       .print dc
c         32       .print tran
c         33       .print ac
c         34       .print noise
c         35       .print distortion
c         36       .plot dc
c         37       .plot tr
c         38       .plot ac
c         39       .plot noise
c         40       .plot distortion
c         41       outputs for dc
c         42       outputs for transient
c         43       outputs for ac
c         44       outputs for noise
c         45       outputs for distortion
c      46-50       <unused>
c
      integer xxor
      dimension lnod(50),lval(50)
      data lnod /10,14,16, 8,15,16,15,16,13, 8,
     1           18,38,27,35, 8, 8,35, 5, 5, 5,
     2            5, 5, 5, 5, 0, 0, 0, 0, 0, 0,
     3           21,21,21,21,21,21,21,21,21,21,
     4            8, 8, 8, 8, 8, 0, 0, 0, 0, 0 /
      data lval / 5, 4, 4, 2, 1, 1, 1, 1, 4, 4,
     1            3, 4, 4,16, 1, 1, 9, 2, 1, 1,
     2           19,55,17,46, 0, 0, 0, 0, 0, 0,
     3            1, 1, 1, 1, 1,17,17,17,17,17,
     4            1, 1, 1, 1, 1, 0, 0, 0, 0, 0 /
      data ndefin /2h.u/
c
c
      anam=aname
      call sizmem(ielmnt,isize)
      locn=ielmnt+isize+2
      if (nsbckt.eq.0) go to 10
      loct=nodplc(isbckt+nsbckt)
      loc=nodplc(loct+3)
      if (loc.ne.0) go to 20
      nodplc(loct+3)=locn
      go to 60
   10 loc=locate(id)
      if (loc.ne.0) go to 20
      locate(id)=locn
      go to 50
c
c  search list for a name match
c
   20 locv=nodplc(loc+1)
      if (xxor(anam,value(locv)).ne.0) go to 30
      if (numalt.ne.0) go to 30
      if (nsbckt.eq.0) go to 25
      if (nodplc(loc-1).ne.id) go to 30
   25 if (nodplc(loc+2).eq.ndefin) go to 200
      if (iforce.eq.0) go to 200
      write (iofile,26) anam
   26 format('0*error*:  above line attempts to redefine ',a8/)
      nogo=1
   30 if (nodplc(loc).eq.0) go to 40
      loc=nodplc(loc)
      go to 20
c
c  reserve space for this element
c
   40 nodplc(loc)=locn
      if (nsbckt.ne.0) go to 60
   50 if (numalt.eq.0) jelcnt(id)=jelcnt(id)+1
   60 loc=locn
      itemp=loc+lnod(id)*nwd4-1
      locv=nxtevn(itemp-1)+1
      itemp=locv-itemp
      ktmp=lnod(id)*nwd4+lval(id)*nwd8+itemp
      call extmem(ielmnt,ktmp)
      locv=(locv-1)/nwd8+1
      iptr=0
      if (nsbckt.eq.0) go to 80
      iptr=id
   80 if (id.le.24) nodplc(loc+lnod(id)-2)=numalt
      nodplc(loc-1)=iptr
      nodplc(loc)=0
      nodplc(loc+1)=locv
      value(locv)=anam
c
c  background storage
c
  100 nodplc(loc+2)=ndefin
      nword=lnod(id)-4
      if (id.le.24) nword=nword-1
      if (nword.lt.1) go to 120
      call zero4(nodplc(loc+3),nword)
  120 nword=lval(id)-1
      if (nword.lt.1) go to 200
      call zero8(value(locv+1),nword)
      if ((id.ge.21).and.(id.le.24)) call undefi(value(locv+1),nword)
c
c  exit
c
  200 return
      end
