#----------------------------------------------------------------
# Revised Makefile (19/02/2014 PS)
# How to call:
# calling as normal
#     > make
# will trigger the default compiler options which include optimisation.
#
# You can now also call using
#     > make CO=debug
# for extra warnings, gprof and gdb output, exception trapping
# at runtime, and bounds-checking.
# This option is slow (about 2x slower than make all)
#
#     > make CO=debug2, debug3, pedantic
# offer further levels of checks in case of problems
#
#     > make new
# simply calls clean then all to force a re-build.
#
#     > (sudo) make install
# places the files in the standard UNIX directories.
#
# I have also included similar options for ifort. Since I have
# the compiler here, and it is potentially significantly faster,
# it may be useful when it comes time to do science.
#
#----------------------------------------------------------------

FC=gfortran
LD=gfortran

# set prefix depending on OS
OS := $(shell uname)
ifeq ($(OS),Darwin)
  PREFIX=/usr/local
else
  PREFIX=/usr
endif

# get the version from the debian changelog if this is a package, and from the git log otherwise
VERSION := $(shell if [ -e debian/ ]; then dpkg-parsechangelog -S version; elif [ "`command -v git`" != "" ]; then git describe --always --tags --dirty; else echo "2.4"; fi)

FFLAGS+=-cpp -DPREFIX=\"$(PREFIX)\" -DVERSION=\"$(VERSION)\"
LDFLAGS+=
MANDIR=$(DESTDIR)$(PREFIX)/share/man/man1

ifeq ($(MESSAGES),yes)
    FFLAGS += -DCO=\"$(CO)\"
endif

ifeq ($(FC),gfortran)
    FFLAGS += -ffree-line-length-0 -Jsource/ -fopenmp
  ifeq ($(CO),debug)
    FFLAGS += -fbounds-check -Wall -Wuninitialized
  else ifeq ($(CO),debug2)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized
  else ifeq ($(CO),debug3)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -ffpe-trap=zero,overflow,invalid,underflow,denormal -fbacktrace -fcheck=all
  else ifeq ($(CO),pedantic)
    FFLAGS += -g -pg -fbounds-check -Wall -Wuninitialized -Werror -pedantic -ffpe-trap=zero,overflow,invalid,underflow,denormal
  else
    FFLAGS += -O3 -fno-backtrace
  endif
endif

ifeq ($(FC),ifort)
  FFLAGS += -module source/ -openmp
  LD=ifort
  ifeq ($(CO),debug)
    FFLAGS += -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback
  else ifeq ($(CO),pedantic)
    FFLAGS += -pg -g -check bounds -check uninit -warn all -warn nodeclarations -WB -zero -traceback -std
  else
    FFLAGS += -fast
  endif
endif

.PHONY: all clean install uninstall new htmlmanual

all: neat

new:
	$(MAKE) clean
	$(MAKE) all

%.o: %.f90
	$(FC) $(FFLAGS) $< -c -o $@

neat: source/globals.o source/types.o source/oii_diagnostics.o source/hydrogen.o source/extinction.o source/recombination_lines.o source/helium.o source/equib_routines.o source/quicksort.o source/functions.o source/filereading.o source/abundances.o source/commandline.o source/linefinder.o source/weights.o source/output.o source/neat.o
	$(LD) $(LDFLAGS) $(FFLAGS) -o $@ $^ -lcfitsio
	@echo "Installation complete"

source/abundances.o : source/abundances.f90 source/globals.o source/hydrogen.o source/oii_diagnostics.o source/types.o source/types.o source/extinction.o source/recombination_lines.o source/helium.o source/filereading.o source/equib_routines.o source/types.o
source/commandline.o : source/commandline.f90 source/functions.o source/globals.o
source/equib_routines.o : source/equib_routines.f90 source/helium.o source/types.o source/globals.o
source/extinction.o : source/extinction.f90 source/hydrogen.o source/globals.o source/types.o
source/filereading.o : source/filereading.f90 source/functions.o source/globals.o source/types.o source/types.o
source/functions.o : source/functions.f90 source/quicksort.o source/types.o source/types.o
source/globals.o : source/globals.f90
source/helium.o : source/helium.f90 source/types.o source/globals.o
source/hydrogen.o : source/hydrogen.f90 source/globals.o source/types.o
source/linefinder.o : source/linefinder.f90 source/quicksort.o source/types.o source/functions.o source/globals.o
source/neat.o : source/neat.f90 source/output.o source/commandline.o source/functions.o source/weights.o source/oii_diagnostics.o source/helium.o source/hydrogen.o source/linefinder.o source/recombination_lines.o source/filereading.o source/types.o source/filereading.o source/quicksort.o source/extinction.o source/types.o source/globals.o source/types.o
source/oii_diagnostics.o : source/oii_diagnostics.f90 source/globals.o
source/output.o : source/output.f90 source/globals.o source/functions.o
source/quicksort.o : source/quicksort.f90 source/globals.o
source/recombination_lines.o : source/recombination_lines.f90 source/globals.o
source/types.o : source/types.f90 source/globals.o
source/weights.o : source/weights.f90 source/globals.o source/filereading.o source/types.o

clean:
	rm -f neat source/*.o source/*.mod man/neat.html

install: neat
	test -e $(DESTDIR)$(PREFIX)/share/neat || mkdir -p $(DESTDIR)$(PREFIX)/share/neat
	test -e $(DESTDIR)$(PREFIX)/share/doc/neat/examples || mkdir -p $(DESTDIR)$(PREFIX)/share/doc/neat/examples
	test -e $(DESTDIR)$(PREFIX)/bin || mkdir -p $(DESTDIR)$(PREFIX)/bin
	test -e $(MANDIR) || mkdir -p $(MANDIR)
	install -m 644 Atomic-data/*.* $(DESTDIR)$(PREFIX)/share/neat
	install -m 644 utilities/complete_line_list $(DESTDIR)$(PREFIX)/share/neat
	install -m 644 utilities/plot.sh $(DESTDIR)$(PREFIX)/share/neat
	install -m 644 config/default.cfg $(DESTDIR)$(PREFIX)/share/neat
	install -m 644 examples/*.dat $(DESTDIR)$(PREFIX)/share/doc/neat/examples
	install neat $(DESTDIR)$(PREFIX)/bin
	install -m 644 man/neat.1 $(MANDIR)
	test -e $(DESTDIR)$(PREFIX)/share/bash-completion/completions || mkdir -p $(DESTDIR)$(PREFIX)/share/bash-completion/completions
	install -m 644 source/bashcompletion $(DESTDIR)$(PREFIX)/share/bash-completion/completions/neat
	gzip -f $(MANDIR)/neat.1
	@echo "Installation complete"

uninstall:
	rm -rf $(DESTDIR)$(PREFIX)/share/neat $(DESTDIR)$(PREFIX)/share/doc/neat
	rm -f $(DESTDIR)$(PREFIX)/bin/neat
	rm -f $(DESTDIR)$(PREFIX)/share/bash-completion/completions/neat
	rm -f $(MANDIR)/neat.1.gz

htmlmanual:
	groff -m mandoc -Thtml man/neat.1 > man/neat.html
