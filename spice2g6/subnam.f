      subroutine subnam(loce)
      implicit double precision (a-h,o-z)
c
c     this routine constructs the names of elements added as a result of
c subcircuit expansion.  the full element names are of the form
c                  name.xn. --- xd.xc.xb.xa
c where 'name' is the nominal element name, and the 'x'*s denote the
c sequence of subcircuit calls (from top or circuit level down through
c nested subcircuit calls) which caused the particular element to be
c added.  at present, spice restricts all element names to be 8 charac-
c ters or less.  therefore, the name used consists of the leftmost 8
c characters of the full element name, with the rightmost character
c replaced by an asterisk ('*') if the full element name is longer than
c 8 characters.
c
c spice version 2g.6  sccsid=blank 3/15/83
      common /blank/ value(200000)
      integer nodplc(64)
      complex cvalue(32)
      equivalence (value(1),nodplc(1),cvalue(1))
c
c
      data ablank, aper, astk / 1h , 1h., 1h* /
c
c  construct subcircuit element name
c
      if (nodplc(loce-1).eq.0) go to 100
      locve=nodplc(loce+1)
      loc=loce
      nchar=0
      sname=ablank
      achar=ablank
   10 locv=nodplc(loc+1)
      elname=value(locv)
      do 20 ichar=1,8
      call move(achar,1,elname,ichar,1)
      if (achar.eq.ablank) go to 30
      if (nchar.eq.8) go to 40
      nchar=nchar+1
      call move(sname,nchar,achar,1,1)
   20 continue
   30 loc=nodplc(loc-1)
      if (loc.eq.0) go to 60
      if (nchar.eq.8) go to 40
      nchar=nchar+1
      call move(sname,nchar,aper,1,1)
      go to 10
c
c  name is longer than 8 characters:  flag with asterisk
c
   40 call move(sname,8,astk,1,1)
   60 value(locve)=sname
c
c  finished
c
  100 return
      end
