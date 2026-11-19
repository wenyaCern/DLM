.PHONY: all estimate smoke run check
GENERATOR := ../../other/build/kernel_generation/build_adaptive_pipi_coulomb_kernel
CHECKER := ../../other/build/kernel_generation/check_adaptive_kernel
KERNEL_OUT ?= /scratch8/Wenya/pipi_coulomb_q24_r72_nonuniform.root
SMOKE_OUT := ../../other/build/kernel_generation/smoke_farm_kernel.root

all:
	$(MAKE) -C ../../kernel_generation all

estimate: all
	$(GENERATOR) $(KERNEL_OUT) --farm --estimate

smoke: all
	$(GENERATOR) $(SMOKE_OUT) --farm --smoke
	$(CHECKER) $(SMOKE_OUT)

run: all
	@test -d $(dir $(KERNEL_OUT)) || { echo "Missing output directory: $(dir $(KERNEL_OUT))"; exit 1; }
	@test ! -e $(KERNEL_OUT) || { echo "Refusing to overwrite: $(KERNEL_OUT)"; exit 1; }
	$(GENERATOR) $(KERNEL_OUT) --farm

check: all
	$(CHECKER) $(KERNEL_OUT) --metadata
