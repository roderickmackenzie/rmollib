hpcinc=-I/apps/gsl/1.8/include/
hpclib=-L/apps/gsl/1.8/lib/
link=-L/apps/gsl/1.8/lib/
#-llzo2
debug=-ggdb
# -pg -g
opp=-O0 $(hpcinc)
main:
	gcc pdbtest.c rmolib.a -o pdbtest -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)
	gcc pdbgenbox.c rmolib.a -o pdbgenbox -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)
	gcc pdbevap.c rmolib.a -o pdbevap -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)
	gcc pdbrandomize.c rmolib.a -o pdbrandomize -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)
	gcc pdbtruncate.c rmolib.a -o pdbtruncate -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)
	gcc pdbcount.c rmolib.a -o pdbcount -lm -lgsl  -lgslcblas $(debug) $(opp) -Wall $(hpclib) $(link)




