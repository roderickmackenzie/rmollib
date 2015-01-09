hpcinc=-I/apps/gsl/1.8/include/ -I/lib/
hpclib=-L/apps/gsl/1.8/lib/
link= -lm -L/apps/gsl/1.8/lib/
#-llzo2
debug=-ggdb
# -pg -g
opp=-O0 $(hpcinc)
main: vec.o util.o save.o random.o progress.o load.o molops.o system.o  evap.o expbuf.o  lzo.o
	ar rvs rmolib.a vec.o util.o save.o random.o progress.o load.o molops.o system.o  evap.o expbuf.o lzo.o
	mv rmolib.a ../

lzo.o:  minilzo.c
	gcc -c  minilzo.c -o lzo.o $(debug) $(opp) -Wall

expbuf.o: expbuf.c
	gcc -c expbuf.c -o expbuf.o $(debug) $(opp) -Wall


evap.o: evap.c
	gcc -c evap.c -o evap.o $(debug) $(opp) -Wall

molops.o: molops.c
	gcc -c molops.c -o molops.o $(debug) -O2 -Wall

system.o: system.c
	gcc -c system.c -o system.o $(debug) $(opp) -Wall

vec.o:	vec.c
	gcc -c vec.c -o vec.o $(debug) $(opp) -Wall


util.o:	util.c
	gcc -c util.c -o util.o $(debug) $(opp)

save.o:	save.c
	gcc -c save.c -o save.o $(debug) $(opp) -Wall

random.o: random.c
	gcc -c random.c -o random.o $(debug) $(opp) -Wall


progress.o: progress.c
	gcc -c progress.c -o progress.o $(debug) $(opp)

load.o: load.c
	gcc -c load.c -o load.o $(debug) $(opp) -Wall


clean:
	rm *.o
