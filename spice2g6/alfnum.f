      subroutine alfnum(number,string,ipos)
      implicit double precision (a-h,o-z)
c
c     this routine converts number into character form, storing the
c characters in the character array string, beginning with the position
c indicated by ipos.
c
c **** note that the 'ipos' variable is changed to indicate the position
c      of the next unwritten character.  this could clobber constants if
c      ipos is not a variable in the calling program
c
      dimension string(1)
      dimension adigit(10)
      data adigit / 1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9 /
      data aminus / 1h- /
c
c
      num=number
c
c  check for number < 0
c
      if (num.ge.0) go to 10
      num=-num
c...  negative number:  insert minus sign
      call move(string,ipos,aminus,1,1)
      ipos=ipos+1
c
c  convert number one digit at a time, in reverse order
c
   10 istart=ipos
   20 numtmp=num/10
      idigit=num-numtmp*10
      call move(string,ipos,adigit(idigit+1),1,1)
      ipos=ipos+1
      num=numtmp
      if (num.ne.0) go to 20
      istop=ipos-1
c
c  now reverse the order of the digits
c
   30 if (istop.le.istart) go to 40
      call move(tmpdgt,1,string,istart,1)
      call move(string,istart,string,istop,1)
      call move(string,istop,tmpdgt,1,1)
      istart=istart+1
      istop=istop-1
      go to 30
c
c  conversion complete
c
   40 return
      end
