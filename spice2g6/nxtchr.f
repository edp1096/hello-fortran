      integer function nxtchr(int)
      implicit double precision (a-h,o-z)
c
c     this routine advances the current line scan pointer one column
c     and checks whether or not the next character is a delimiter
c
c spice version 2g.6  sccsid=line 3/15/83
      common /line/ achar,afield(15),oldlin(15),kntrc,kntlim
c
      dimension adelim(5)
      data adelim / 1h , 1h,, 1h=, 1h(, 1h) /
      data ablnk / 1h  /
      data ichar /0/
c
c  advance scan pointer (kntrc)
c
      kntrc=kntrc+1
      if (kntrc.gt.kntlim) go to 30
      call move(achar,1,afield,kntrc,1)
    5 do 10 i=1,5
      if (achar.eq.adelim(i)) go to 20
   10 continue
c
c  non-delimiter
c
      nxtchr=1
      return
c
c  delimiter
c
   20 nxtchr=0
      return
c
c  end-of-line
c
   30 nxtchr=-1
      achar=ablnk
      return
      end
