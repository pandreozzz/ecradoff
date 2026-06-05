# make can be invoked using "make PROFILE=<prof>" in which case your
# local configuration parameters will be obtained from
# Makefile_include.<prof>
ifndef PROFILE
$(info *** No "PROFILE" variable provided, assuming "gfortran")
PROFILE = gfortran
endif

ifndef ECRADOFF_DP
$(info *** Default ECRADOFF_DP=0 for single-precision fields in ecradoff)
ECRADOFF_DP = 0
endif

# Include a platform-specific makefile that defines FC, FCFLAGS and
# LIBS
include	ecrad/Makefile_include.$(PROFILE)

SOURCES = interp2d.f90 interpvert.f90 interplut.f90

INCLUDE_DIR = locals/lib/ #libs/shared

OBJECTS := $(addprefix locals/lib/,$(SOURCES:.f90=.so))
# $(addprefix libs/shared/,$(SOURCES:.f90=.so))

all: $(OBJECTS)

FCSHAREDFLAGS = -fpic -fopenmp -shared


# If CPPFLAGS is empty and PROFILE is intel or intel_atos, set it to fpp, else fail
ifeq ($(CPPFLAGS),)
	ifeq ($(PROFILE),intel)
		CPPFLAGS = -fpp
	else ifeq ($(PROFILE),intel_atos)
		CPPFLAGS = -fpp
	else
		$(error "CPPFLAGS empty!")
	endif
endif

ifeq ($(ECRADOFF_DP),1)
	CPPFLAGS += -DECRADOFF_DP
endif

#libs/shared/%.so: libs/src/%.f90
locals/lib/%.so: src/libs/%.f90
	$(FC) $(FCFLAGS) $(FCSHAREDFLAGS) $(CPPFLAGS) -o $@ $<

clean:
	rm -f *.o $(OBJECTS)
