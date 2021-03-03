source= mod_calculate_rates.f90 mod_global_functions.f90 mod_global_variables.f90 dvode_f90_m.f90 mod_read_model.f90 mod_read_rate06.f90 chem_rate06_dvode.f90 mod_save_results.f90 mod_run_dvode.f90


# Format:
#formFLAG = -noauto -Wall

# real underflow/ integer overflow:
#flowFLAG = -check underflow -check overflow

# Optimierung:
#optFLAG = -Ofast  -march=native -ffree-line-length-512 

# Memory debugging
debugFLAG = -g -pg -Wall -fbacktrace -Wunderflow -Woverflow -ffree-line-length-512 

FLAGS = $(debugFLAG) $(optFLAG) $(formFLAG) $(flowFLAG)

monaco: $(source)
	@gfortran $(FLAGS) $(source) -o $@
	@echo make complete

clean:
	rm -rf *~ *.o 
	rm -rf ab
	rm -rf csv
	rm results*
	rm analytics*
	rm *.mod
	rm fort*
