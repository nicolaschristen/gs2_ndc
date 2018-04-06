
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_gs2_optimisation

#Define this tests level
MY_TEST_LEVEL=2
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	$(call SKIPCOMMAND,test_level is ${TEST_LEVEL}: this test is level ${MY_TEST_LEVEL})
else
all: test_gs2_optimisation 
	$(call RUNCOMMAND,test_gs2_optimisation,test_gs2_optimisation.in)
	$(call RUNCOMMAND,../../../gs2,test_gs2_optimisation_run.in)
endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_gs2_optimisation.o:  ${TEST_MODS} 

test_gs2_optimisation: test_gs2_optimisation.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

