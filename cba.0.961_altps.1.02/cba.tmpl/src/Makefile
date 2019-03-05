########################################################################
# NOTE: replace 'PROJECT' and 'project' below
#       according to your project name
########################################################################

########################################################################
# NOTE: need to set environmental variables 'CBA_HOME'
#       e.g.,
#         export CBA_HOME=$HOME/cba
########################################################################

CPP	= g++
CPPFLAGS = 
AR = ar
ARFLAGS = -rv
.SUFFIXES: .o .cpp .h
.cpp.o:
	$(CPP) $(CPPFLAGS) $(INCLUDE) -c $<

LIBDIR = ../lib
INCDIR = ../include
BINDIR = ../bin

CBA_LIBDIR = $(CBA_HOME)/lib
CBA_INCDIR = $(CBA_HOME)/include
CBA_BINDIR = $(CBA_HOME)/bin

LIBCBA = $(CBA_LIBDIR)/libcba.a
LIBMAIN = $(CBA_LIBDIR)/libcbamain.a
LIBPROJECT = $(LIBDIR)/libcbaproject.a

INCLUDE = -I$(CBA_INCDIR)
LIBS = -lm -L$(LIBDIR) -lcbaproject\
	-L$(CBA_LIBDIR) -lcbamain -lcba
LDFLAGS =

# NOTE:
# add utility classes' header and object file names to HEADERS and OBJS;
# add applications' executable file names to APPS;
# e.g.,
# HEADERS = project.h
# OBJS = project.o
# APPS = app_project
HEADERS = 
OBJS = 
APPS = 

all: lib $(APPS)
lib: $(LIBPROJECT)

######################################################################
# applications
######################################################################

#
# NOTE: sample
#app_project: app_project.o $(OBJS) $(LIBCBA) $(LIBMAIN)
#	$(CPP) -o $@ app_project.o $(LIBS)
#app_project.o: app_project.cpp app_project.h $(HEADERS) $(LIBCBA) $(LIBMAIN)

######################################################################
# utilities
######################################################################

$(LIBPROJECT): $(OBJS)
	$(AR) $(ARFLAGS) $@ $?
	-cp -p $(HEADERS) $(INCDIR)
# NOTE: sample
# project.o: project.cpp project.h $(LIBCBA)

install:
	-cp -p $(APPS) $(BINDIR)

clean:
	rm -f *.o *.exe $(APPS)
