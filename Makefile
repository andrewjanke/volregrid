PROGS = volregrid
HEADERS = arb_path_io.h
OBJS = $(PROGS:=.o) $(HEADERS:.h=.o) lex.o

LEX=flex
CC=cc

OPTIONS = -g3 -fullwarn -O3
INCLUDES = -I/usr/local/include
CFLAGS = $(OPTIONS) $(INCLUDES) `gsl-config --cflags`

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lvolume_io -lminc -lnetcdf -lm
LDOPTS = $(LDINCLUDES) $(LDLIBS) `gsl-config --libs`


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OBJS) -o $@ $(OPTIONS) $(LDOPTS)

clean:
	rm -f *.o *~ $(PROGS)
