EXEC   = sfh_test

OBJS   = ./sfh_test.o \
	 ./star_formation_history.o


INCL   = allvars.h proto.h Makefile

CC       = gcc -g -Wall
OPTIMIZE =   
LIBS     =

$(EXEC): $(OBJS) 
	$(CC) $(OPTIMIZE) $(OBJS) $(LIBS)   -o  $(EXEC)  

$(OBJS): $(INCL)

clean:
	rm -f $(OBJS)

tidy:
	rm -f $(OBJS) .$(EXEC)
