include Makefile.SRCS


###################################################

F90 = mpiifort
DEBUG_FLAGS = -init=snan,arrays -check all -check noarg_temp_created -debug-parameters all \
              -traceback -ftrapuv -g -fpe0 -O0
F90FLAGS = -qopenmp -O3  -r8 -g -ftz

LLIB = -lnetcdff
LDFLAGS =  $(F90FLAGS) $(LLIB)

PROG = uEMEPv2

NILUDIR = NILU

%.o: $(NILUDIR)/%.for
	$(F90) $(F90FLAGS) -c $< -o $@

%.o: %.f90
	$(F90) $(F90FLAGS) -c $< -o $@

all:  $(PROG)

$(PROG): dependencies

include dependencies

$(PROG): $(FOBJ)
	$(F90) -o $@ $(FOBJ) $(LDFLAGS)

all:  $(PROG)

clean: 
	rm -f $(PROG) $(FOBJ)
