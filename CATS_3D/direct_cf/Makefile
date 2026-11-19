PROJECT := $(abspath ..)
BUILD := $(PROJECT)/other/build
CMAKE_ARGS ?=

.PHONY: all configure
all: run_CATS3D

configure:
	cmake -S $(PROJECT)/other/build_config -B $(BUILD) $(CMAKE_ARGS)

run_CATS3D: configure
	cmake --build $(BUILD) --target run_CATS3D

# Direct-CATS comparisons that do not fit Data can be built separately.
comparisons: configure
	cmake --build $(BUILD) --target source_compare_CATS3D r20_check_CATS3D calculate_cats3d_r40_q5_octant
