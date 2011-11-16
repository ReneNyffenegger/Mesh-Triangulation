# Fairly generic makefile

# compiler
CC=g++
# compiler flags
CFLAGS=-c -O3
# linker flags
LDFLAGS=
# source files
SRCDIR=.
SRCFILES=meshGeneration.cpp
# executable name
EXECUTABLE=mesh

# test sources
TSTDIR=.
TSTFILES=

####################################################
## EVERYTHING AFTER HERE IS HANDLED AUTOMATICALLY ##
####################################################

SOURCES=$(addprefix $(SRCDIR)/,$(SRCFILES))

TST=$(addprefix $(TSTDIR)/,$(TSTFILES)) $(SOURCES:$(SRCDIR)/main.cpp=)

# test executable
TSTEXECUTABLE=$(addprefix t,$(EXECUTABLE))

# object files (handled automatically)
OBJECTS=$(SOURCES:.cpp=.o)
TSTOBJECTS=$(TST:.cpp=.o)


##############################################################################
all: $(SOURCES) $(EXECUTABLE)

test: $(SOURCES) $(TST) $(TSTEXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $(EXECUTABLE)

$(TSTEXECUTABLE): $(TSTOBJECTS)
	$(CC) $(LDFLAGS) $(TSTOBJECTS) -o $(TSTEXECUTABLE)

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -rf $(OBJECTS) $(EXECUTABLE) $(TSTOBJECTS) $(TSTEXECUTABLE)

