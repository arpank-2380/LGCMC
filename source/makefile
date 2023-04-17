COMPILER = ifort
OUTPUT := ./lgcmc.x
FLAGS = -g  -check all -O0 -traceback 
SOURCE = derived_data_module logo_module variable cubic_eos_module read_input_files site_finder_module calc_ads_en calc_lat_int energy_calculation_module ran_num_gen gcmc_module
EXECUTABLE = main

OBJS = source/$(SOURCE:%=%.o)
ALL_OBJS = $(SOURCE:%=%.o) $(EXECUTABLE:%=%.o)

.PHONY : all
all: $(ALL_OBJS) $(EXECUTABLE)

.PHONY : clean 
clean:
	@echo "Cleaning old objective files and .mod files"
	@\rm -rf *.o *.mod  $(OUTPUT)
%.o: ./%.f90
	@echo "Compiling: " $*
	$(COMPILER) -c $(FLAGS) $< -o $@
        
%: %.o $(ALL_OBJS)
	@echo "   Linking program:   " $a
	$(COMPILER) $(FLAGS) $(ALL_OBJS) -o $(OUTPUT)
