PROGS=regrid_view
OBJS=$(PROGS:=.o) ../arb_path_io.o ../lex.o trackball.o

LEX=flex

OPTIONS = -g3 -O3
INCLUDES = -I/usr/local/mni/include -I. -I../
CFLAGS = $(OPTIONS) $(INCLUDES)

LDINCLUDES = -L/usr/local/mni/lib -L/usr/X11R6/lib
LDLIBS = -lm -lminc -lnetcdf -lglut -lGL -lGLU -lpthread -lX11 -lXext -lXmu
LDOPTS = $(LDINCLUDES) $(LDLIBS)


all: $(PROGS) 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

$(PROGS): $(OBJS)
	$(CC) $(OPTIONS) $(OBJS) -o $@ $(LDOPTS)
	
clean:
	rm -f *.o *~ $(PROGS)

