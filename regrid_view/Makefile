PROGS=regrid_view
OBJS=$(PROGS:=.o) ../arb_path_io.o ../lex.o trackball.o

LEX=flex
CC=cc

OPTIONS = -g3 -fullwarn -O3
INCLUDES = -I/usr/local/include -I. -I../
CFLAGS = $(OPTIONS) $(INCLUDES)

LDINCLUDES = -L/usr/local/lib32
LDLIBS = -lm -lminc -lnetcdf -lglut -lGL -lGLU -lX11 -lXext -lXmu
LDOPTS = $(LDINCLUDES) $(LDLIBS)


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $@ $(LDOPTS)
	
clean:
	rm -f *.o *~ $(PROGS)

