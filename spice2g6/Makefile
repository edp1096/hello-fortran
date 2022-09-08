OPT	= -O3
FFLAGS	= $(OPT) -w
CFLAGS	= $(OPT) -ansi -pedantic
LDFLAGS	= $(OPT) -lm -lf2c
TIME	= time
SIZE	= size
PROG = spice

.f.c:
	f2c -A -c $*.f

OBJ =	spice.c tmpupd.c find.c alter.c title.c dcdcmp.c dcsol.c \
	dmpmat.c indxx.c swapij.c reserv.c matloc.c codgen.c setmem.c \
	locf.c getm4.c getm8.c getm16.c getmx.c avlm4.c avlm8.c \
	avlm16.c relmem.c extmem.c sizmem.c clrmem.c ptrmem.c crunch.c \
	errmem.c memadj.c slpmem.c nxtevn.c nxtmem.c comprs.c memptr.c \
	dmpmem.c memory.c magphs.c xxor.c outnam.c alfnum.c undefi.c \
	getcje.c cpytb4.c cpytb8.c second.c readin.c alias.c keysrc.c \
	extnam.c runcon.c outdef.c card.c getlin.c ushift.c nxtchr.c \
	errchk.c shlsrt.c putnod.c getnod.c subckt.c fndnam.c newnod.c \
	addelt.c lnkref.c subnam.c elprnt.c modchk.c topchk.c setup.c \
	matptr.c reordr.c dctran.c pheadr.c comcof.c trunc.c terr.c \
	sorupd.c sorstp.c iter8.c load.c nlcsrc.c update.c evpoly.c \
	evterm.c nxtpwr.c intgr8.c pnjlim.c diode.c bjt.c fetlim.c \
	limvds.c jfet.c mosfet.c moseq1.c moseq2.c moseq3.c moscap.c \
	cmeyer.c mosq2.c mqspof.c mosq3.c dcop.c sstf.c sencal.c \
	asol.c acan.c cdiv.c cmult.c acdcmp.c acsol.c acload.c noise.c \
	acasol.c dinit.c disto.c ovtpvt.c ntrpl8.c setprn.c setplt.c \
	plot.c scale.c fouran.c unix.c

$(PROG): $(OBJ)

clean:
	/bin/rm -f $(OBJ)
