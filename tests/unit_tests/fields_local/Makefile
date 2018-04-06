
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_fields_local test_fields_local_tuning

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_fields_local test_fields_local_tuning 
	$(call RUNCOMMAND,test_fields_local,test_fields_local.in)
	$(call RUNCOMMAND,test_fields_local_tuning,test_fields_local_tuning.in)

endif

#$(info "TEST_MODS=$(TEST_MODS)")

test_fields_local.o:  ${TEST_MODS} 

test_fields_local: test_fields_local.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

test_fields_local_tuning.o:  ${TEST_MODS}

test_fields_local_tuning: test_fields_local_tuning.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

