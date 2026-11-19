PROJECT := $(abspath ..)
BUILD := $(PROJECT)/other/build
CMAKE_ARGS ?=

.PHONY: help configure gaussian resonance gaussian_all resonance_all
help:
	@echo "Choose one target: gaussian, resonance, gaussian_all, resonance_all"

configure:
	cmake -S $(PROJECT)/other/build_config -B $(BUILD) $(CMAKE_ARGS)

gaussian: configure
	cmake --build $(BUILD) --target fit_CATS3D_derived_full

resonance: configure
	cmake --build $(BUILD) --target fit_resonance6291_mt0

gaussian_all: configure
	cmake --build $(BUILD) --target fit_CATS3D fit_CATS3D_oo fit_CATS3D_derived_full

resonance_all: configure
	cmake --build $(BUILD) --target fit_CATS3D_resonance fit_resonance6291_mt0
