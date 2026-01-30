# ==============================================================================
# Professional Makefile for Coupled/Russian Simulation
# Supports SERIAL and PARALLEL (MUMPS + OpenMP) builds
# ==============================================================================

# ==============================================================================
# Configuration Variables
# ==============================================================================

# Detect if 'parallel' is in the command line goals
ifneq (,$(filter parallel,$(MAKECMDGOALS)))
    PARALLEL := 1
endif

# Parallelization mode (can also be set via PARALLEL=1 on command line)
PARALLEL      ?= 0

# Compiler settings
CXX           := g++
CXXSTANDARD   := -std=c++17

# Build configuration
BUILD_TYPE    ?= release
VERBOSE_MODE  ?= 0

# ==============================================================================
# GetFEM Installation Path Detection
# ==============================================================================

ifeq ($(PARALLEL),1)
    # Parallel version: use getfem-parallel installation
    GETFEM_PREFIX ?= $(HOME)/getfem/getfem-5.4.4-parallel
    DEFINES       := -DUSE_MUMPS
else
    # Serial version: use standard getfem installation
    GETFEM_PREFIX ?= $(HOME)/getfem/getfem-5.4.4-standard
    DEFINES       :=
endif

# Try to auto-detect if not manually set
# (User can override with: make GETFEM_PREFIX=/custom/path)
ifndef GETFEM_PREFIX_OVERRIDE
    # Check if the default path exists, otherwise try to find it
    ifeq ($(wildcard $(GETFEM_PREFIX)),)
        # Try common locations
        ifneq ($(wildcard /usr/local/include/getfem),)
            GETFEM_PREFIX := /usr/local
        else ifneq ($(wildcard /usr/include/getfem),)
            GETFEM_PREFIX := /usr
        else ifneq ($(wildcard $(HOME)/getfem),)
            # Use first found in home directory
            GETFEM_PREFIX := $(shell find $(HOME)/getfem -maxdepth 1 -type d -name "getfem-*" | head -n1)
        endif
    endif
endif

GETFEM_INC    := $(GETFEM_PREFIX)/include
GETFEM_LIB    := $(GETFEM_PREFIX)/lib

# ==============================================================================
# Project Directories
# ==============================================================================

SRC_DIR       := src
INC_DIR       := include
OBJ_DIR       := obj

# ==============================================================================
# Source Files
# ==============================================================================

ALL_SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
SOURCES       := $(filter-out $(SRC_DIR)/main_%.cc,$(ALL_SRC_FILES))
OBJECTS       := $(SOURCES:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

# ==============================================================================
# Main File Detection
# ==============================================================================

# Look for main_coupled
MAIN_COUPLED_SRC :=
MAIN_COUPLED_OBJ :=
MAIN_COUPLED_DIR :=
ifneq ($(wildcard main_coupled.cc),)
    MAIN_COUPLED_SRC := main_coupled.cc
    MAIN_COUPLED_OBJ := $(OBJ_DIR)/main_coupled.o
    MAIN_COUPLED_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_coupled.cc),)
    MAIN_COUPLED_SRC := main_coupled.cc
    MAIN_COUPLED_OBJ := $(OBJ_DIR)/main_coupled.o
    MAIN_COUPLED_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_coupled.cpp),)
    MAIN_COUPLED_SRC := main_coupled.cpp
    MAIN_COUPLED_OBJ := $(OBJ_DIR)/main_coupled.o
    MAIN_COUPLED_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_coupled.cpp),)
    MAIN_COUPLED_SRC := main_coupled.cpp
    MAIN_COUPLED_OBJ := $(OBJ_DIR)/main_coupled.o
    MAIN_COUPLED_DIR := $(SRC_DIR)
endif

# Look for main_russian_doll or main_russian
MAIN_RUSSIAN_SRC :=
MAIN_RUSSIAN_OBJ :=
MAIN_RUSSIAN_DIR :=
ifneq ($(wildcard main_russian_doll.cc),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cc
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian_doll.o
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian_doll.cc),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cc
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian_doll.o
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian_doll.cpp),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cpp
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian_doll.o
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian_doll.cpp),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cpp
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian_doll.o
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian.cc),)
    MAIN_RUSSIAN_SRC := main_russian.cc
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian.o
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian.cc),)
    MAIN_RUSSIAN_SRC := main_russian.cc
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian.o
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian.cpp),)
    MAIN_RUSSIAN_SRC := main_russian.cpp
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian.o
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian.cpp),)
    MAIN_RUSSIAN_SRC := main_russian.cpp
    MAIN_RUSSIAN_OBJ := $(OBJ_DIR)/main_russian.o
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
endif

# Extract base names for targets
MAIN_COUPLED := $(basename $(notdir $(MAIN_COUPLED_SRC)))
MAIN_RUSSIAN := $(basename $(notdir $(MAIN_RUSSIAN_SRC)))

# Executables
TARGET_COUPLED := $(MAIN_COUPLED)
TARGET_RUSSIAN := $(MAIN_RUSSIAN)

# Set default target
ifneq ($(MAIN_COUPLED_SRC),)
    DEFAULT_TARGET := coupled
else ifneq ($(MAIN_RUSSIAN_SRC),)
    DEFAULT_TARGET := russian
else
    DEFAULT_TARGET := 
endif

# ==============================================================================
# Compiler Flags
# ==============================================================================

# Include paths
INCLUDES      := -I$(INC_DIR) -I$(GETFEM_INC)

# Add MPI includes when in parallel mode
ifeq ($(PARALLEL),1)
    INCLUDES  += -I/usr/lib/x86_64-linux-gnu/openmpi/include \
                 -I/usr/lib/x86_64-linux-gnu/openmpi/include/openmpi
endif

# Base compiler flags
CXXFLAGS_BASE := $(CXXSTANDARD) $(INCLUDES) -Wall -Wextra

# Add verbose flag if enabled
ifeq ($(VERBOSE_MODE),1)
    DEFINES += -DVERBOSE
endif

# Build-type specific flags
ifeq ($(BUILD_TYPE),debug)
    CXXFLAGS_OPT := -g -O0 -DDEBUG
else ifeq ($(BUILD_TYPE),release)
    CXXFLAGS_OPT := -O3 -DNDEBUG -march=native
else ifeq ($(BUILD_TYPE),profile)
    CXXFLAGS_OPT := -g -O2 -pg
else
    $(error Invalid BUILD_TYPE: $(BUILD_TYPE). Use 'debug', 'release', or 'profile')
endif

# Warning suppression flags
WARN_SUPPRESS := -Wno-comment \
                 -Wno-unused-but-set-variable \
                 -Wno-unused-parameter \
                 -Wno-pragmas \
                 -Wno-unknown-pragmas \
                 -Wno-deprecated-declarations \
                 -Wno-misleading-indentation \
                 -Wno-cast-function-type

# Combined compiler flags
CXXFLAGS      := $(CXXFLAGS_BASE) $(CXXFLAGS_OPT) $(DEFINES) $(WARN_SUPPRESS)

# ==============================================================================
# Linker Flags and Libraries
# ==============================================================================

# Linker flags
LDFLAGS       := -L$(GETFEM_LIB)

# Add MPI library path when in parallel mode
ifeq ($(PARALLEL),1)
    LDFLAGS   += -L/usr/lib/x86_64-linux-gnu/openmpi/lib
endif

# Libraries - order matters!
LDLIBS        := -lgetfem

# CRITICAL: If using parallel GetFEM, -lgomp MUST come next
# Auto-detect if GetFEM path contains "parallel" OR if PARALLEL=1
ifneq (,$(findstring parallel,$(GETFEM_PREFIX)))
    LDLIBS    += -lgomp
else ifeq ($(PARALLEL),1)
    LDLIBS    += -lgomp
endif

# MUMPS libraries (only in parallel mode)
ifeq ($(PARALLEL),1)
    LDLIBS    += -ldmumps \
                 -lzmumps \
                 -lcmumps \
                 -lsmumps \
                 -lmumps_common \
                 -lpord \
                 -lscotch \
                 -lscotcherr \
                 -lmetis
endif

# SuperLU (sparse direct solver - always needed)
LDLIBS        += -lsuperlu

# Qhull (computational geometry)
LDLIBS        += -lqhull_r

# BLAS/LAPACK - try pkg-config first, then manual detection
LAPACK_LIBS   := $(shell pkg-config --libs lapack 2>/dev/null)
BLAS_LIBS     := $(shell pkg-config --libs blas 2>/dev/null)

ifeq ($(LAPACK_LIBS),)
    # Detect OpenBLAS or use standard libraries
    ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblas.so*),)
        BLAS_LIBS := -lopenblas
    else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libopenblas.so*),)
        BLAS_LIBS := -lopenblas
    else
        BLAS_LIBS := -lblas
    endif
    LAPACK_LIBS := -llapack
endif

LDLIBS        += $(LAPACK_LIBS) $(BLAS_LIBS)

# Fortran runtime (needed for MUMPS/LAPACK)
LDLIBS        += -lgfortran

# MPI C++ libraries (CRITICAL - must come after all other libraries)
ifeq ($(PARALLEL),1)
    LDLIBS    += -lmpi_cxx -lmpi -lopen-rte -lopen-pal
endif

# System libraries
LDLIBS        += -lpthread -lm

# ==============================================================================
# Phony Targets
# ==============================================================================

.PHONY: all coupled russian clean distclean help info directories \
        check-getfem list-mains parallel serial

# Default target
all: $(DEFAULT_TARGET)

# ==============================================================================
# Convenience Targets
# ==============================================================================

# Serial and parallel are marker targets
# The actual parallelism is controlled by MAKECMDGOALS detection above

# If only 'serial' is specified, build default target
ifeq ($(MAKECMDGOALS),serial)
serial: all
else
serial:
endif
	@# Serial marker

# If only 'parallel' is specified, build default target
ifeq ($(MAKECMDGOALS),parallel)
parallel: all
else
parallel:
endif
	@# Parallel marker

# ==============================================================================
# Validation Targets
# ==============================================================================

check-getfem:
	@echo "==> Checking GetFEM installation..."
	@if [ ! -d "$(GETFEM_INC)" ]; then \
		echo "ERROR: GetFEM include directory not found: $(GETFEM_INC)"; \
		echo "Please set GETFEM_PREFIX to your GetFEM installation:"; \
		echo "  make GETFEM_PREFIX=/path/to/getfem"; \
		exit 1; \
	fi
	@if [ ! -d "$(GETFEM_LIB)" ]; then \
		echo "ERROR: GetFEM library directory not found: $(GETFEM_LIB)"; \
		exit 1; \
	fi
	@echo "✓ GetFEM found at: $(GETFEM_PREFIX)"
ifeq ($(PARALLEL),1)
	@echo "✓ Parallel mode enabled (MUMPS + OpenMP)"
else
	@echo "✓ Serial mode"
endif

list-mains:
	@echo "==> Detected main source files:"
ifneq ($(MAIN_COUPLED_SRC),)
	@echo "  Coupled: $(MAIN_COUPLED_DIR)/$(MAIN_COUPLED_SRC)"
endif
ifneq ($(MAIN_RUSSIAN_SRC),)
	@echo "  Russian: $(MAIN_RUSSIAN_DIR)/$(MAIN_RUSSIAN_SRC)"
endif
ifeq ($(MAIN_COUPLED_SRC)$(MAIN_RUSSIAN_SRC),)
	@echo "  No main files found! Looking for:"
	@echo "    - main_coupled.{cc,cpp} in . or src/"
	@echo "    - main_russian_doll.{cc,cpp} in . or src/"
	@echo "    - main_russian.{cc,cpp} in . or src/"
endif
	@echo ""
	@echo "==> Available targets:"
ifneq ($(MAIN_COUPLED_SRC),)
	@echo "  make coupled   -> builds $(TARGET_COUPLED)"
endif
ifneq ($(MAIN_RUSSIAN_SRC),)
	@echo "  make russian   -> builds $(TARGET_RUSSIAN)"
endif

# ==============================================================================
# Main Build Targets
# ==============================================================================

coupled: check-getfem directories $(TARGET_COUPLED)
	@echo ""
	@echo "==> Build complete: $(TARGET_COUPLED)"
ifeq ($(PARALLEL),1)
	@echo "==> Parallel mode (MUMPS + OpenMP)"
	@echo "==> Run with: export OMP_NUM_THREADS=4; ./$(TARGET_COUPLED) -f data_file.txt"
else
	@echo "==> Serial mode"
	@echo "==> Run with: ./$(TARGET_COUPLED) -f data_file.txt"
endif
	@echo ""

russian: check-getfem directories $(TARGET_RUSSIAN)
	@echo ""
	@echo "==> Build complete: $(TARGET_RUSSIAN)"
ifeq ($(PARALLEL),1)
	@echo "==> Parallel mode (MUMPS + OpenMP)"
	@echo "==> Run with: export OMP_NUM_THREADS=4; ./$(TARGET_RUSSIAN) -f data_file.txt"
else
	@echo "==> Serial mode"
	@echo "==> Run with: ./$(TARGET_RUSSIAN) -f data_file.txt"
endif
	@echo ""

# ==============================================================================
# Linking Rules
# ==============================================================================

$(TARGET_COUPLED): $(OBJECTS) $(MAIN_COUPLED_OBJ)
	@echo "==> Linking $@"
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(TARGET_RUSSIAN): $(OBJECTS) $(MAIN_RUSSIAN_OBJ)
	@echo "==> Linking $@"
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# ==============================================================================
# Compilation Rules
# ==============================================================================

# Compile .cc files from src directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cc files from root directory (for main files)
$(OBJ_DIR)/%.o: %.cc | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cpp files from src directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Compile .cpp files from root directory (for main files)
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@

# Create necessary directories
directories: $(OBJ_DIR)

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# ==============================================================================
# Utility Targets
# ==============================================================================

clean:
	@echo "==> Cleaning object files and binaries"
	@rm -rf $(OBJ_DIR)
	@rm -f $(TARGET_COUPLED) $(TARGET_RUSSIAN)

distclean: clean
	@echo "==> Removing all generated files"
	@find . -name "*~" -delete
	@find . -name "*.o" -delete
	@find . -name "core" -delete
	@rm -f *.mm *.vtk
	@rm -rf output_vtk

info:
	@echo "===================================================================="
	@echo "Build Configuration"
	@echo "===================================================================="
	@echo "Compiler:          $(CXX)"
	@echo "C++ Standard:      $(CXXSTANDARD)"
	@echo "Build Type:        $(BUILD_TYPE)"
	@echo "Verbose Mode:      $(VERBOSE_MODE)"
	@echo "Parallel Mode:     $(PARALLEL)"
ifeq ($(PARALLEL),1)
	@echo "  - MUMPS solver:  ENABLED"
	@echo "  - OpenMP:        ENABLED (via GetFEM)"
	@echo "  - MPI library:   ENABLED (for MUMPS)"
else
	@echo "  - Serial build"
endif
	@echo "--------------------------------------------------------------------"
	@echo "GetFEM Prefix:     $(GETFEM_PREFIX)"
	@echo "GetFEM Include:    $(GETFEM_INC)"
	@echo "GetFEM Library:    $(GETFEM_LIB)"
	@echo "--------------------------------------------------------------------"
ifneq ($(MAIN_COUPLED_SRC),)
	@echo "Main coupled:      $(MAIN_COUPLED_DIR)/$(MAIN_COUPLED_SRC)"
	@echo "  -> Output:       ./$(TARGET_COUPLED)"
else
	@echo "Main coupled:      NOT FOUND"
endif
ifneq ($(MAIN_RUSSIAN_SRC),)
	@echo "Main russian:      $(MAIN_RUSSIAN_DIR)/$(MAIN_RUSSIAN_SRC)"
	@echo "  -> Output:       ./$(TARGET_RUSSIAN)"
else
	@echo "Main russian:      NOT FOUND"
endif
	@echo "Default target:    $(DEFAULT_TARGET)"
	@echo "--------------------------------------------------------------------"
	@echo "Source files:      $(words $(SOURCES)) files"
	@echo "Object files:      $(words $(OBJECTS)) objects"
	@echo "--------------------------------------------------------------------"
	@echo "CXXFLAGS:          $(CXXFLAGS)"
	@echo "LDFLAGS:           $(LDFLAGS)"
	@echo "LDLIBS:            $(LDLIBS)"
	@echo "===================================================================="

help:
	@echo "===================================================================="
	@echo "Makefile for Coupled/Russian Poroelasticity Simulation"
	@echo "===================================================================="
	@echo ""
	@echo "QUICK START:"
	@echo "  make                      # Serial build (default)"
	@echo "  make parallel             # Parallel build (default target)"
	@echo "  make parallel russian     # Parallel russian build"
	@echo "  make PARALLEL=1 russian   # Alternative parallel syntax"
	@echo "  make russian              # Serial russian build"
	@echo ""
	@echo "===================================================================="
	@echo "Available Targets"
	@echo "===================================================================="
	@echo "  make [all]           - Build default target ($(DEFAULT_TARGET))"
ifneq ($(MAIN_COUPLED_SRC),)
	@echo "  make coupled         - Build $(TARGET_COUPLED)"
endif
ifneq ($(MAIN_RUSSIAN_SRC),)
	@echo "  make russian         - Build $(TARGET_RUSSIAN)"
endif
	@echo "  make serial          - Explicit serial build"
	@echo "  make parallel        - Parallel build (MUMPS + OpenMP)"
	@echo "  make clean           - Remove object files and binaries"
	@echo "  make distclean       - Remove all generated files"
	@echo "  make info            - Display build configuration"
	@echo "  make help            - Display this help message"
	@echo "  make check-getfem    - Verify GetFEM installation"
	@echo "  make list-mains      - List detected main source files"
	@echo ""
	@echo "===================================================================="
	@echo "Build Options (set via command line)"
	@echo "===================================================================="
	@echo "  PARALLEL=<0|1>       - Enable parallel mode (MUMPS + OpenMP)"
	@echo "                         0 = serial (default), 1 = parallel"
	@echo "  BUILD_TYPE=<type>    - Set build type"
	@echo "                         debug, release (default), profile"
	@echo "  VERBOSE_MODE=<0|1>   - Enable VERBOSE preprocessor flag"
	@echo "  GETFEM_PREFIX=<path> - Set GetFEM installation prefix"
	@echo "                         Auto-detected by default"
	@echo ""
	@echo "===================================================================="
	@echo "Usage Examples"
	@echo "===================================================================="
	@echo "  # Serial builds"
	@echo "  make                              # Default serial"
	@echo "  make coupled                      # Coupled simulation"
	@echo "  make russian                      # Russian doll simulation"
	@echo ""
	@echo "  # Parallel builds (RECOMMENDED for large problems)"
	@echo "  make parallel                     # Parallel coupled"
	@echo "  make parallel russian             # Parallel russian"
	@echo "  make PARALLEL=1 coupled -j4       # Parallel compile with 4 cores"
	@echo ""
	@echo "  # Debug builds"
	@echo "  make BUILD_TYPE=debug             # Serial debug"
	@echo "  make parallel BUILD_TYPE=debug    # Parallel debug"
	@echo ""
	@echo "  # Custom GetFEM location"
	@echo "  make GETFEM_PREFIX=/custom/path"
	@echo ""
	@echo "  # Verbose output"
	@echo "  make VERBOSE_MODE=1"
	@echo ""
	@echo "  # Clean output (suppress Boost pragma warnings):"
	@echo "  chmod +x clean-make"
	@echo "  ./clean-make parallel russian -j4"
	@echo ""
	@echo "===================================================================="
	@echo "Running the Simulation"
	@echo "===================================================================="
	@echo "  # Serial mode:"
	@echo "  ./main_coupled -f data_file.txt"
	@echo ""
	@echo "  # Parallel mode (set thread count):"
	@echo "  export OMP_NUM_THREADS=4"
	@echo "  export OPENBLAS_NUM_THREADS=4"
	@echo "  ./main_coupled -f data_file.txt"
	@echo ""
	@echo "===================================================================="
	@echo "Required Libraries (Parallel Mode)"
	@echo "===================================================================="
	@echo "  - GetFEM (compiled with --enable-openmp --enable-mumps)"
	@echo "  - MUMPS (parallel sparse direct solver)"
	@echo "  - OpenMPI (MPI library)"
	@echo "  - OpenBLAS or BLAS/LAPACK (linear algebra)"
	@echo "  - SuperLU (sparse solver fallback)"
	@echo "  - SCOTCH, METIS (graph partitioning)"
	@echo "  - Qhull (computational geometry)"
	@echo ""
	@echo "  Install on Ubuntu/Debian:"
	@echo "    sudo apt-get install libopenmpi-dev libmumps-dev"
	@echo "    sudo apt-get install libopenblas-dev liblapack-dev"
	@echo "    sudo apt-get install libscotch-dev libmetis-dev"
	@echo "    sudo apt-get install libsuperlu-dev libqhull-dev"
	@echo "===================================================================="

# ==============================================================================
# Dependency Generation
# ==============================================================================

DEPS := $(OBJECTS:.o=.d)
ifneq ($(MAIN_COUPLED_OBJ),)
    DEPS += $(MAIN_COUPLED_OBJ:.o=.d)
endif
ifneq ($(MAIN_RUSSIAN_OBJ),)
    DEPS += $(MAIN_RUSSIAN_OBJ:.o=.d)
endif

-include $(DEPS)

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: %.cc | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: %.cpp | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@