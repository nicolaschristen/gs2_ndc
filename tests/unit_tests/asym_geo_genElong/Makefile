
sinclude ../Makefile.tests_common

CLEANCOMMAND+= test_asym_geo_genElong

#Define this tests level
MY_TEST_LEVEL=1
VAL:= $(shell expr ${MY_TEST_LEVEL} \<= ${TEST_LEVEL})
ifeq "$(VAL)" "0"
all:
	@echo "==> Skipping test in ${PWD}/$(TEST_FOLD) (Test_level ${TEST_LEVEL}, This test ${MY_TEST_LEVEL})."
else
all: test_asym_geo_genElong
	$(call RUNCOMMAND,test_asym_geo_genElong,)

endif

ogyropsi.dat EQDSK:
	tar -xjf data.tar.bz2

test_asym_geo_genElong.o:  ${TEST_MODS} 

test_asym_geo_genElong: test_asym_geo_genElong.o $(TEST_MODS)
	$(LD) $(LDFLAGS) -o $@ $^ $(LIBS)

