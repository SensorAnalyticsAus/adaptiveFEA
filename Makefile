PROGRAM  = femgs1 
PROGRAM2 = adpgs1
PROGRAM3 = faopgs2
PROGRAM4 = plotp

# Set CC compilation flags

CC = gcc
CFLAGS = -O    
LDLIBS = -lm

# Set FC compilation

FC = gfortran

OBJS = fem.o renum.o jcc.o 

OBJS2 = adp.o

OBJS3 = aom1.o aow1.o

HEADERS3 = incl.c inclx.c


$(PROGRAM) : $(PROGRAM2) $(OBJS)
	$(CC) $(OBJS) $(LDLIBS) -o $(PROGRAM) 

$(PROGRAM2) : $(OBJS2)
	$(CC) $(OBJS2) $(LDLIBS) -o $(PROGRAM2)

$(PROGRAM3) : $(OBJS3)
	$(CC) $(OBJS3) $(LDLIBS) -o $(PROGRAM3) 

$(OBJS3) : $(HEADERS3) 

$(PROGRAM4):
	$(FC) plotp01.f -o $(PROGRAM4)

all: $(PROGRAM) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4)

.PHONY: clean distclean
clean : 
	rm -f core *.o *.*% *.str *.me *_.dat *.fem *.err *.adp *.inf
distclean:
	rm -f core *.o *.*% *.str *.me *_.dat *.fem *.err *.adp *.inf
	rm -f adpgs1 femgs1 faopgs2 plotp 
