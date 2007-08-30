
# Get the machine type and kernel name. Piping the output of uname
# through (tr '[A-Z]' '[a-z]') converts the string to lower case

ARCH := $(shell (uname -m | tr '[A-Z]' '[a-z]'))
OS := $(shell (uname -s | tr '[A-Z]' '[a-z]'))
ARCHOS := $(ARCH)_$(OS)

ifeq ($(ARCHOS),i686_linux)
  CC = gcc
  CFLAGS =  -Wall -Wunused-variable -ansi
  DEBUGFLAGS = -g -DDEBUG
  OPTFLAGS = -O                 # -pg for PROFILING
  DEFINES = -DLINUX -DHASHTABLE 
  LDFLAGS := -lm
endif

# METIS DIRECTORY

METISDIR = $(HOME)/packages/metis-4.0
METISINC = $(METISDIR)/Lib


# STANDARD DIRECTORY

TOPDIR = $(HOME)
DEPINCS := -I$(METISINC)

# DEVELOP DIRECTORY

TOPDIR = $(HOME)/develop
DEPINCS := -I$(METISINC)

# ---------------- DO NOT EDIT BELOW THIS LINE --------------------------------

.PHONY: all
all: debug opt


LIBDIR = ./lib

INCDIR = -I./include $(DEPINCS) 

ifeq ($(PAR),1)
srcdirs := src/util src/base src/hilev src/misc src/par
else
srcdirs := src/util src/base src/hilev src/misc
endif
incdirs := include

# Sources with path

srcs1 := $(foreach dir,$(srcdirs),$(wildcard $(dir)/*.c))


# Headers with path

hdrs := $(foreach dir,$(incdirs),$(wildcard $(dir)/*.h))


# Sources without path

srcs := $(notdir $(srcs1))


# Object files without path

objs := $(srcs:.c=.o)

# Object file directories (-d is for debug)

objdir-d := obj/$(ARCHOS)-d
objdir   := obj/$(ARCHOS)

# Debug object files with path

obj2-d := $(addprefix $(objdir-d)/,$(objs))


# Optimized object files with path

obj2 := $(addprefix $(objdir)/,$(objs))


# Directory for libraries

libdir := lib/$(ARCHOS)


VPATH = $(objdir-d):$(objdir):$(dirs)




# DEBUG VERSION OF THE LIBRARY

debug: $(objdir-d) $(libdir) $(obj2-d)
	 ar rcv $(LIBDIR)/$(ARCHOS)/libmstk-d.a $(obj2-d)

$(objdir-d):
	mkdir -p $(objdir-d)

obj/$(ARCHOS)-d/%.o: src/util/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/par/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@



# OPTIMIZED VERSION OF THE LIBRARY

opt : $(objdir) $(libdir) $(obj2)
	ar rcv $(LIBDIR)/$(ARCHOS)/libmstk.a $(obj2)

$(objdir):
	mkdir -p $(objdir)

obj/$(ARCHOS)/%.o: src/util/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

obj/$(ARCHOS)/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

obj/$(ARCHOS)/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


obj/$(ARCHOS)/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


obj/$(ARCHOS)/%.o: src/par/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


$(libdir):
	mkdir -p $(libdir)


# CLEAN

clean :
	rm -f $(obj2-d) $(obj2) $(LIBDIR)/$(ARCHOS)/*.a


# DEPEND

depend:
	makedepend  $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.
