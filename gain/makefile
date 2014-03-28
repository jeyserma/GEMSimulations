
GARFIELD =/afs/cern.ch/user/j/jeyserma/garfieldpp

OBJDIR = $(GARFIELD)/Object
SRCDIR = $(GARFIELD)/Source
INCDIR = $(GARFIELD)/Include
HEEDDIR = $(GARFIELD)/Heed
LIBDIR = $(GARFIELD)/Library


# Compiler flags
CFLAGS = -Wall -Wextra -Wno-long-long \
`root-config --cflags` \
-O3 -fno-common -c \
-I$(INCDIR) -I$(HEEDDIR)

# Debug flags
#CFLAGS += -g

LDFLAGS = `root-config --glibs` -lgfortran -lm
LDFLAGS += -L$(LIBDIR) -lGarfield
#LDFLAGS += -g

gain: gain.c
	$(CXX) $(CFLAGS) gain.c
	$(CXX) -o gain gain.o $(LDFLAGS)
	rm gain.o
