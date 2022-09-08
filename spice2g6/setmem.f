      subroutine setmem(ipntr,ksize)
      implicit double precision (a-h,o-z)
c
c     this routine performs dynamic memory management.  it is used in
c     spice2, and useable in any program.
c
c     memory is managed within an array selected by the calling program.
c     one may either dimension this array to the 'maxmem' size, or more
c     desirably, find the address of the first available word of memory
c     above your program, and dimension your array to '1'.  passing the
c     address of the first data word available permits the manager to
c     use 'illegal' indices into the data area.
c
c     this routine must have access to an integer function called 'locf'
c     which returns the address of its argument.  addresses as used by this
c     program refer to 'integer' addresses, not byte addresses.
c
c entry points:
c      setmem - set initial memory
c      getm4  - get block for table of integers
c      getm8  - get block for table of floating point variables
c      getm16 - get block for table of complex variables
c      relmem - release part of block
c      extmem - extend size of existing block
c      sizmem - determine size of existing block
c      clrmem - release block
c      ptrmem - reset memory pointer
c      crunch - force memory compaction
c      avlm4  - amount of space available (integers)
c      avlm8  - amount of space available (real)
c      avlm16 - amount of space available (complex)
c
c calling sequences:
c      call setmem(imem(1),maxmem)
c      call setmem(imem(1),maxmem,kfamwa)  cdc machines running under
c                                          calidoscope kfamwa is the address
c                                          of the first available word
c      call getm4 (ipntr,blksiz)  where blksize is the number of entries
c      call getm8 (ipntr,blksiz)
c      call getm16(ipntr,blksiz)
c      call relmem(ipntr,relsiz)
c      call extmem(ipntr,extsiz)  extsiz is the number of entries to be added
c      call sizmem(ipntr,blksiz)
c      call clrmem(ipntr)
c      call ptrmem(ipntr1,ipntr2)
c      call avlm4(ispace)
c      call avlm8(ispace)
c      call avlm16(ispace)
c      call crunch
c      call slpmem(ipntr,slpsiz)  express desire for *slpsiz* extra entries
c
c
c general comments:
c      for each block which is allocated, a multi-word entry is maintained
c in a table kept in high memory, of the form
c
c        word      contents
c        ----      --------
c
c          1       index of imem(.) into origin of block
c                    i.e. contents of pointer (used for error check)
c          2       block size (in words)
c          3       number of words in use
c          4       address of variable containing block origin
c          5       number of words used per table entry
c          6       slop size (in words)
c
c      all allocated blocks are an 'even' (nxtevn) number of words in length,
c where a 'word' is the storage unit required for an 'integer' variable.
c      since block repositioning may be necessary, the convention that
c only one variable contain a block origin should be observed.
c      for *getmem*, *ipntr* is set such that *array(ipntr+1)* is the
c first word of the allocated block.  'ipntr' is set to address the first
c entry of the table when used with the appropriate variable type, i.e.,
c nodplc(ipntr+1), value(ipntr+1), or cvalue(ipntr+1).
c      for *clrmem*, *ipntr* is set to 'invalid' to enable rapid detection
c of an attempt to use a cleared block.
c      if any fatal errors are found, a message is printed and a flag
c set inhibiting further action until *setmem* is called.  (in this
c context, insufficient memory is considered a fatal error.)
c      throughout this routine, *ldval* always contains the subscript of
c the last addressable word of memory, *memavl* always contains the
c number of available words of memory, *numblk* always contains the
c number of allocated blocks, and istack(*loctab* +1) always contains
c the first word of the block table.
c
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
c spice version 2g.6  sccsid=memmgr 3/15/83
      common /memmgr/ cpyknt,istack(1),lorg,icore,maxcor,maxuse,memavl,
     1   ldval,numblk,loctab,ltab,ifwa,nwoff,ntab,maxmem,memerr,nwd4,
     2   nwd8,nwd16
      dimension ipntr(1)
c
      logical memptr
      complex cvalue(32)
      equivalence (value(1),cvalue(1))
      external locf
c
c...  approximate time required to copy *nwords* integer values
c
c  nxtevn rounds the number up to the next 'even' value.  the value
c  used for this 'even' number is the smallest number into which one
c  can divide nwd4,nwd8,and nwd16.
c
c
c  nxtmem  returns next higher memory size
c
c
c
c***  setmem - set initial memory
c
      nwd4=1
      nwd8=locf(value(2))-locf(value(1))
      nwd16=locf(cvalue(2))-locf(cvalue(1))
      memerr=0
      nevn=nxtevn(1)
c     check that nxtevn function returns a number divisible by
c     nwd4, nwd8, nwd16; also check that the memory increment
c     nxtmem(.) is an integer multiple of nxtevn(1)
      icheck=mod(nevn,nwd4)+mod(nevn,nwd8)+mod(nevn,nwd16)+
     1  mod(nxtmem(1),nevn)
      if(icheck.eq.0) go to 2
      memerr=1
      call errmem(6,memerr,ipntr(1))
    2 cpyknt=0.0d0
      ifamwa=locf(ipntr(1))
      maxmem=ksize
      ntab=nxtevn(6)
c... add 'lorg' to an address and you get the 'istack' index to that word
      lorg=1-locf(istack(1))
      ifwa=ifamwa+lorg-1
      nwoff=locf(ipntr(1))+lorg-1
      icore=nxtmem(1)
c... don't take chances, back off from 'end of memory' by nxtevn(1)
      ldval=ifwa+nxtmem(1)-nxtevn(1)
      memavl=ldval-ntab-ifwa
      maxcor=0
      maxuse=0
      call memory
      if(memerr.ne.0) call errmem(6,memerr,ipntr(1))
      numblk=1
      loctab=ldval-ntab
      istack(loctab+1)=0
      istack(loctab+2)=memavl
      istack(loctab+3)=0
      istack(loctab+4)=-1
      istack(loctab+5)=1
      istack(loctab+6)=0
      return
      end
