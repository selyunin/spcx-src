EXTENSION=.lin64
C_LIBS=/home/kaergel/amv/trunk/c
PREFIX=aaf_lib
LIB=libaaf.lin64.a
LIBNAME   = $(LIB:%.a=%)
SHAREDLIB = $(LIBNAME).so
DEBUG = 

MAKE = make -j
MAKEFILE = Makefile

CCC = g++
CCFLAGS += -DLIN64 -Wall -DUSE_AAF_EXTENSIONS -Wno-deprecated -O3
CCFLAGS += -Weffc++ $(DEBUG)

AR = ar
ARFLAGS = ruv
RANLIB = ranlib

RM = rm -f 

DOXYGEN = doxygen

SOURCES.c =	aa_rounding.c 	aa_interval.c 	aa_aaftrigo.c 	aa_aafapprox.c \
		aa_aafarithm.c 	aa_aafcommon.c 	aa_exceptions.c	aa_aafspecial.c

SOURCES.h =	aa.h aa_aaf.h 	aa_interval.h 	aa_rounding.h	aa_exceptions.h 

SOURCES = $(SOURCES.h)  $(SOURCES.c)

OBJECTS = $(SOURCES.c:%.c=%$(EXTENSION).o)

OBJECTS_SO = $(SOURCES.c:%.c=%$(EXTENSION).so.o)

default: all

all:	 $(LIB) $(SHAREDLIB)

static:  $(LIB)

shared:  $(SHAREDLIB)

debug:
	@$(MAKE) -f $(MAKEFILE) cleaner
	@$(MAKE) -f $(MAKEFILE) DEBUG='-g'

$(LIB): $(SOURCES) $(OBJECTS)
	@echo "\n - Compiling $(LIB) -\n"
	@rm -f $(LIB)
	$(AR) $(ARFLAGS) $(LIB) $(OBJECTS)
	$(RANLIB) $(LIB)

%$(EXTENSION).o: %.c $(SOURCES.h)
	$(CCC) $(CCFLAGS) -c $(@:%$(EXTENSION).o=%.c) -o $@                                             

$(SHAREDLIB): $(OBJECTS_SO)
	gcc -shared -Wall -o $@  $(OBJECTS_SO)

%$(EXTENSION).so.o: %.c $(SOURCES.h)
	$(CCC) $(CCFLAGS) -fPIC -c $(@:%$(EXTENSION).so.o=%.c) -o $@                                             

clean::
	$(RM) $(OBJECTS) $(OBJECTS_SO)

cleaner:: clean
	$(RM) $(LIB) $(SHAREDLIB) config.status
	$(RM) -r doxygen

doxygen: $(SOURCES) doxygen.*
	$(DOXYGEN) doxygen.cfg
	@echo "doxygen done."

