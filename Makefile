PROGS = volregrid
HEADERS = arb_path_io.h minc_support.h
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o) lex.o

LEX=flex
CC=gcc

# WARN = -Wall
WARN = -Wall -Wtraditional -Wshadow -Wpointer-arith -Wcast-qual \
       -Wcast-align -Wconversion -Waggregate-return -Wstrict-prototypes \
       -Wmissing-prototypes -Wmissing-declarations -Wredundant-decls \
       -Winline -pedantic -Wunused -Wunused-parameter

OPTIONS = $(WARN) -O3
INCLUDES = -I/usr/local/mni/include
CFLAGS = $(OPTIONS) $(INCLUDES) `gsl-config --cflags`

LDINCLUDES = -L/usr/local/mni/lib -L/usr/local/lib
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS) `gsl-config --libs`


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
