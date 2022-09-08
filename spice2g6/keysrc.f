      subroutine keysrc(keytab,lentab,tstwrd,index)
      implicit double precision (a-h,o-z)
      double precision keytab
c
c     this routine searches the keyword table 'keytab' for the possible
c entry 'tstwrd'.  abbreviations are considered as matches.
c
      dimension keytab(lentab)
      integer xxor
      data ablnk / 1h  /
c
c
      index=0
      lenwrd=0
      achar=ablnk
      do 10 i=1,8
      call move(achar,8,tstwrd,i,1)
      if (achar.eq.ablnk) go to 20
      lenwrd=lenwrd+1
   10 continue
c
   20 if (lenwrd.eq.0) go to 40
      tstchr=ablnk
      call move(tstchr,8,tstwrd,1,1)
   30 index=index+1
      if (index.gt.lentab) go to 40
      akey=ablnk
      call move(akey,1,keytab(index),1,lenwrd)
      if (xxor(akey,tstwrd).eq.0) go to 50
      go to 30
c
   40 index=-1
   50 return
      end
