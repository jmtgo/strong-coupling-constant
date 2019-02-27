

all: f2ave.exe testvar.exe

files = source/f2ave.o source/getcovar.o source/readf2.o source/w3cehiq200.o\
 source/f2swim.o source/h12000pdf.o source/toblockdiag.o source/w3cphiq200.o\
 source/fillarrays.o source/w2cehiq200.o source/wlcehiq200.o\
 source/fitf2.o source/outf2.o source/w2cphiq200.o source/wlcphiq200.o \
 source/f2rnd.o source/xscorrcme.o \
 source/f2nchiq200.o  source/f3nchiq200.o   source/flnchiq200.o \
 source/dsdQ2.o source/statrecalc.o source/simf2.o source/readrea.o source/readtable.o \
 source/recovar.o source/chi2res.o  source/flxf3proc.o  source/writef2fl.o \
 source/charm.o source/asymp.o \
 source/cg0f2lh.o source/cg1f2lh.o \
 source/cq1f2lh.o source/daind.o source/dilog.o \
 source/hblg_ca.o source/hblq_h.o source/hbtg_ca.o \
 source/hbtq_h.o source/hlg_ca.o source/hlg_cf.o \
 source/hlq_h.o source/htg_ca.o source/htg_cf.o \
 source/htq_h.o source/locate.o source/thresh.o \
 source/hlq_l.o source/htq_l.o source/f2charm.o\
 source/mrst2004f3.o source/cteq5.o\
 source/redx.o source/zpdfn.o

f2ave.exe: $(files) include/*.inc
	g77  -fno-globals -fno-automatic -finit-local-zero -o f2ave.exe \
	${files} -L/afs/cern.ch/sw/lcg/external/cernlib/2006a/i686-slc5-gcc34-opt/lib -lmathlib -lpacklib -lkernlib -llapack3 -lblas -L /usr/lib/gcc/i386-redhat-linux/4.1.2  -lgfortran

testvar.exe: source/chi2var.o source/testvar.o include/*.inc
	g77 -fno-globals -fno-automatic -finit-local-zero -o testvar.exe \
	source/chi2var.o source/testvar.o

testrecovar.exe: source/readrea.o source/readtable.o source/testrecovar.o source/recovar.o
	g77 -fno-globals -fno-automatic -finit-local-zero -o testrecovar.exe \
	source/readrea.o source/readtable.o source/testrecovar.o source/recovar.o \
	-L/cern/pro/lib -lmathlib -lpacklib -lkernlib -llapack3 -lblas

.f.o:
	g77 -fno-globals    -fno-automatic -finit-local-zero -c  $< ; \
        mv *.o source/


