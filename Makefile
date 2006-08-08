# For now we are just compiling on i686_linux
#ARCHOS := $(shell ($(HOME)/bin/getarch))

ARCHOS = i686_linux

ifeq ($(ARCHOS),i686_linux)
  CC = gcc
  CFLAGS =  -Wall -Wunused-variable -ansi
  DEBUGFLAGS = -g -DDEBUG
  OPTFLAGS = -O                 # -pg for PROFILING
  DEFINES = -DLINUX 
  LDFLAGS := -lm
endif

# METIS DIRECTORY

METISDIR = $(HOME)/packages/metis-4.0
METISINC = $(METISDIR)/Lib


# STANDARD DIRECTORY

TOPDIR = $(HOME)
MSTKUTIL_VER = 1.4
DEPINCS := -I$(TOPDIR)/include/mstkutil-$(MSTKUTIL_VER) -I$(METISINC)

# DEVELOP DIRECTORY

TOPDIR = $(HOME)/develop
MSTKUTIL_VER = 1.5dev
DEPINCS := -I$(TOPDIR)/mstkutil/$(MSTKUTIL_VER)/include -I$(METISINC)

# ---------------- DO NOT EDIT BELOW THIS LINE --------------------------------

.PHONY: all
all: debug opt


LIBDIR = ./lib

INCDIR = -I./include $(DEPINCS) 

ifeq ($(PAR),1)
srcdirs := src/base src/hilev src/misc src/par
else
srcdirs := src/base src/hilev src/misc
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


# Debug object files with path

obj2-d := $(addprefix obj/$(ARCHOS)-d/,$(objs))


# Optimized object files with path

obj2 := $(addprefix obj/$(ARCHOS)/,$(objs))


VPATH = obj/$(ARCHOS)-d:obj/$(ARCHOS):$(dirs)




# DEBUG VERSION OF THE LIBRARY

debug: $(obj2-d)
	 ar rcv $(LIBDIR)/$(ARCHOS)/libmstk-d.a $^

obj/$(ARCHOS)-d/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

obj/$(ARCHOS)-d/%.o: src/par/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@



# OPTIMIZED VERSION OF THE LIBRARY

opt : $(obj2)
	ar rcv $(LIBDIR)/$(ARCHOS)/libmstk.a $^ 

obj/$(ARCHOS)/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

obj/$(ARCHOS)/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


obj/$(ARCHOS)/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@


obj/$(ARCHOS)/%.o: src/par/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@



# CLEAN

clean :
	rm -f $(obj2-d) $(obj2) $(LIBDIR)/$(ARCHOS)/*.a


# DEPEND

depend:
	makedepend  $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.
