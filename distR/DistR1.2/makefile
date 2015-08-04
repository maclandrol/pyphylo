SRCDIR = ./src
BINDIR = ./bin
INCDIR = ./include
PROGNAME=distR

distR:
	export PROGNAME
	cd src; $(MAKE) clean
	cd src; $(MAKE) all
	mkdir -p $(BINDIR)
	cp $(SRCDIR)/$(PROGNAME) $(BINDIR)

clean:
	cd $(SRCDIR) && $(MAKE) clean
