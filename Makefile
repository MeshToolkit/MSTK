
# For now we are just compiling on i686_linux
#ARCHOS := $(shell ($(HOME)/bin/getarch))

ARCHOS = i686_linux

ifeq ($(ARCHOS),i686_linux)
  CC = gcc
  CFLAGS =  -Wall -Wno-unused-variable -ansi
  DEBUGFLAGS = -g -DDEBUG
  OPTFLAGS = -O -pg
  DEFINES = -DLINUX 
  LDFLAGS := -lm
endif

TOPDIR = $(HOME)/develop
MSTK_UTIL := $(TOPDIR)/mstkutil/1.2

# ---------------- DO NOT EDIT BELOW THIS LINE --------------------------------

LIBDIR = ./lib
MSTK_UTIL_INC = -I$(MSTK_UTIL)/include

INCDIR = -I./include $(MSTK_UTIL_INC) 

dirs := src/base src/hilev src/misc
temp := $(foreach dir,$(dirs),$(wildcard $(dir)/*.c))
hdrs := $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
files := $(temp) $(foreach dir,$(dirs),$(wildcard $(dir)/*.h))
srcs := $(notdir $(temp))
objs := $(srcs:.c=.o)
obj2-d := $(addprefix obj/$(ARCHOS)-d/,$(objs))
obj2 := $(addprefix obj/$(ARCHOS)/,$(objs))
VPATH = obj/$(ARCHOS)-d:obj/$(ARCHOS):$(dirs)


# DEBUG VERSION OF THE LIBRARY

mstk : $(obj2-d)
	 ar rcv $(LIBDIR)/$(ARCHOS)/libmstk-d.a $^

obj/$(ARCHOS)-d/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c src/base/$< -o obj/$(ARCHOS)-d/$@

obj/$(ARCHOS)-d/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c src/hilev/$< -o obj/$(ARCHOS)-d/$@

obj/$(ARCHOS)-d/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c $< -o $@

%.o: %.cc
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(DEBUGFLAGS) -c src/misc/$< -o obj/$(ARCHOS)-d/$@

# OPTIMIZED VERSION OF THE LIBRARY

opt : $(obj2)
	ar rcv $(LIBDIR)/$(ARCHOS)/libmstk.a $^ 

obj/$(ARCHOS)/%.o: src/base/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c src/base/$< -o obj/$(ARCHOS)/$@

obj/$(ARCHOS)/%.o: src/hilev/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c src/hilev/$< -o obj/$(ARCHOS)/$@

obj/$(ARCHOS)/%.o: src/misc/%.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(DEFINES) $(INCDIR) $(CFLAGS) $(OPTFLAGS) -c src/misc/$< -o obj/$(ARCHOS)/$@

clean :
	rm $(obj2-d) $(obj2) $(LIBDIR)/$(ARCHOS)/*.a

print: $(files)
	ppc $?
	touch print

printh: $(hdrs)
	ppc $?
	touch  printh

depend:
	makedepend  $(INCDIR) $(temp)

# DO NOT DELETE THIS LINE -- make depend depends on it.
