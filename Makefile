PROGS = volregrid
HEADERS = arb_path_io.h
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o) lex.o

LEX=flex
CC=cc

OPTIONS = -pedantic -Wall -O3
INCLUDES = -I/usr/local/mni/include
CFLAGS = $(OPTIONS) $(INCLUDES) `gsl-config --cflags`

LDINCLUDES = -L/usr/local/mni/lib
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS) `gsl-config --libs`


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
