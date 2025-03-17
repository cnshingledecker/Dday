# Source files
SOURCES = mod_global_variables.f90 mod_global_functions.f90 mod_calculate_rates.f90 dvode_f90_m.f90 \
          mod_read_model.f90 mod_read_rate06.f90 mod_save_results.f90 mod_run_dvode.f90 chem_rate06_dvode.f90

# Object files (derived from source files)
OBJECTS = $(SOURCES:.f90=.o)

# Compiler and flags
#FC = ifort
 FC = gfortran  # Uncomment for gfortran
#FLAGS = -g -pg -O2 -fprotect-parens -fp-model=strict -march=native -module . -diag-disable=10448  # ifort options with module path and warning suppression
 FLAGS = -g -pg -O2 -march=native -ffree-line-length-512 -J .  # gfortran options with module path

# Target executable
TARGET = monaco

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJECTS)
	$(FC) $(FLAGS) -o $@ $(OBJECTS)      # TAB required here
	@echo "Make complete"                # TAB required here

# Compilation rules with dependencies
mod_global_variables.o: mod_global_variables.f90
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_global_functions.o: mod_global_functions.f90 mod_global_variables.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_calculate_rates.o: mod_calculate_rates.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

dvode_f90_m.o: dvode_f90_m.f90
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_read_rate06.o: mod_read_rate06.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_read_model.o: mod_read_model.f90 mod_global_variables.o mod_global_functions.o mod_read_rate06.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_save_results.o: mod_save_results.f90 mod_global_variables.o mod_global_functions.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

mod_run_dvode.o: mod_run_dvode.f90 mod_global_variables.o mod_global_functions.o mod_calculate_rates.o dvode_f90_m.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

chem_rate06_dvode.o: chem_rate06_dvode.f90 mod_global_variables.o mod_global_functions.o mod_calculate_rates.o \
                     dvode_f90_m.o mod_read_model.o mod_read_rate06.o mod_save_results.o mod_run_dvode.o
	$(FC) $(FLAGS) -c $< -o $@           # TAB required here

# Clean up
clean:
	rm -rf *~ *.o *.mod ab csv fort.* results* analytics*  # TAB required here

# Phony targets
.PHONY: all clean
