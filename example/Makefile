# Get the machine type and kernel name. Piping the output of uname
# # through (tr '[A-Z]' '[a-z]') converts the string to lower case


ARCH := $(shell (uname -m | tr '[A-Z]' '[a-z]'))
OS := $(shell (uname -s | tr '[A-Z]' '[a-z]'))
ARCHOS := $(ARCH)_$(OS)


BASEDIR := $(HOME)/develop

ifeq ($(ARCHOS),i686_linux)
CC = gcc
CFLAGS = -g -Wall
OCFLAGS = -O -Wall
DEFINES = -DDEBUG -DDEBUG2
endif
ifeq ($(ARCHOS),x86_64_linux)
CC = gcc
CFLAGS = -g -Wall
OCFLAGS = -O -Wall
DEFINES = -DDEBUG -DDEBUG2
endif


DEPS := mstk/1.7 

INCDIR   := -I./include $(addsuffix /include,$(addprefix -I$(BASEDIR)/,$(DEPS)))
LIBDIR := $(addsuffix /lib/$(ARCHOS),$(addprefix -L$(BASEDIR)/,$(DEPS)))

LIBS := $(LIBDIR) \
	-lmstk-d 

OPTLIBS := $(LIBDIR) \
	-lmstk

OTHERLIBS := -lm

OBJ := obj/$(ARCHOS)-d
OPTOBJ := obj/$(ARCHOS)

dir_code:=src
dir_obj:=obj
temp  := $(foreach dir,$(dir_code),$(wildcard $(dir)/*.c))
srcs  := $(notdir $(temp))
objs  := $(srcs:.c=.o)
obj2  := $(addprefix obj/$(ARCHOS)-d/,$(objs))
optobj2  := $(addprefix obj/$(ARCHOS)/,$(objs))
VPATH=src

INSTALLDIR := ./bin


debug : $(obj2)
	$(CC) $(LDFLAGS) $(INCDIR) $^ $(LDOBJS) $(LIBS) $(OTHERLIBS)  -o $(INSTALLDIR)/$(ARCHOS)/example-d
	@echo
	@echo
	@echo "Executable is $(INSTALLDIR)/$(ARCHOS)/example-d"
	@echo "DEBUGGING IS ON"
	@echo

$(OBJ)/%.o : src/%.c
	$(CC) $(DEFINES) $(INCDIR)  $(CFLAGS) $(EXTRACFLAGS) -c $< -o $@

opt : $(optobj2)
	$(CC) $(LDFLAGS)  $(INCDIR) $^ $(LDOBJS) $(OPTLIBS) $(OTHERLIBS) -o $(INSTALLDIR)/$(ARCHOS)/example 
	@echo
	@echo
	@echo "Executable is $(INSTALLDIR)/$(ARCHOS)/example"
	@echo "OPTIMIZED VERSION, DEBUGGING IS OFF"
	@echo


$(OPTOBJ)/%.o : src/%.c
	$(CC) $(DEFINES) $(INCDIR)  $(OCFLAGS) $(EXTRACFLAGS) -c $< -o $@


clean :
	/bin/rm $(obj2) $(optobj2)

depend:
	makedepend $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.



