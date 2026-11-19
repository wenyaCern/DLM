CXX ?= c++
CATS_INSTALL := $(abspath ../../install)
CATS_CONFIG := $(CATS_INSTALL)/bin/cats-config
ROOT_CFLAGS := $(shell root-config --cflags)
ROOT_LIBS := $(shell root-config --libs)
CATS_INCLUDE := $(shell $(CATS_CONFIG) --incdir)
CATS_LIBDIR := $(shell $(CATS_CONFIG) --libdir)
BIN_DIR := ../other/build/kernel_generation
CPPFLAGS += -I../fit/shared -I../fit/common -I.. -I$(CATS_INCLUDE) $(ROOT_CFLAGS)
CXXFLAGS += -std=c++17 -O2
LDLIBS += -L$(CATS_LIBDIR) -Wl,-rpath,$(CATS_LIBDIR) -lCATS -lMinuit $(ROOT_LIBS)

.PHONY: all smoke clean legacy
all: $(BIN_DIR)/build_adaptive_pipi_coulomb_kernel $(BIN_DIR)/check_adaptive_kernel

$(BIN_DIR):
	mkdir -p $@

$(BIN_DIR)/build_adaptive_pipi_coulomb_kernel: build_adaptive_pipi_coulomb_kernel.cpp | $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< -o $@ $(LDLIBS)

$(BIN_DIR)/check_adaptive_kernel: check_adaptive_kernel.cpp ../fit/shared/AdaptivePipiKernel.cpp ../fit/shared/AdaptivePipiKernel.h | $(BIN_DIR)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) check_adaptive_kernel.cpp ../fit/shared/AdaptivePipiKernel.cpp -o $@ $(LDLIBS)

smoke: all
	$(BIN_DIR)/build_adaptive_pipi_coulomb_kernel $(BIN_DIR)/smoke_adaptive_kernel.root --smoke
	$(BIN_DIR)/check_adaptive_kernel $(BIN_DIR)/smoke_adaptive_kernel.root

# Older kernel builders share the fitter library and remain separate targets.
legacy:
	cmake -S ../other/build_config -B ../other/build
	cmake --build ../other/build --target build_local_source_kernel_CATS3D build_data_kernel_CATS3D convolve_kernel_gauss_CATS3D kernel_cache_check_CATS3D

clean:
	rm -f $(BIN_DIR)/build_adaptive_pipi_coulomb_kernel $(BIN_DIR)/check_adaptive_kernel $(BIN_DIR)/smoke_adaptive_kernel.root
