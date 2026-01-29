# ==============================================================================
# Professional Makefile for Coupled/Russian Simulation
# ==============================================================================

# ==============================================================================
# Configuration Variables
# ==============================================================================

# Compiler settings
USE_MPI       ?= 0
ifeq ($(USE_MPI),1)
    CXX       := mpic++
    DEFINES   := -DUSE_MPI
else
    CXX       := g++
    DEFINES   :=
endif

CXXSTANDARD   := -std=c++17

# Build configuration
BUILD_TYPE    ?= release
VERBOSE_MODE  ?= 0

# External library paths (customize these for your environment)
# First try pkg-config, then fall back to manual path or default
GETFEM_PREFIX = $(HOME)/getfem/getfem-5.4.4-standard
#GETFEM_PREFIX ?= $(shell pkg-config --variable=prefix getfem 2>/dev/null)
#ifeq ($(GETFEM_PREFIX),)
    # Fallback: check common installation directories
#    ifneq ($(wildcard /usr/local/include/getfem),)
   #     GETFEM_PREFIX := /usr/local
 #   else ifneq ($(wildcard /usr/include/getfem),)
  #      GETFEM_PREFIX := /usr
  # else ifneq ($(wildcard $(HOME)/getfem),)
   #     GETFEM_PREFIX := $(HOME)/getfem
 #   else
        # Last resort: use /usr/local and warn user
#        GETFEM_PREFIX := /usr/local
#        $(warning GetFEM not found via pkg-config. Using default: $(GETFEM_PREFIX))
#        $(warning If this is incorrect, set GETFEM_PREFIX manually: make GETFEM_PREFIX=/your/path)
#    endif
#endif

GETFEM_INC    := $(GETFEM_PREFIX)/include
GETFEM_LIB    := $(GETFEM_PREFIX)/lib

# Project directories
SRC_DIR       := src
INC_DIR       := include
OBJ_DIR       := obj

# Source files (only .cc files in src directory, excluding main files)
ALL_SRC_FILES := $(wildcard $(SRC_DIR)/*.cc)
# Filter out any main_* files that might be in src
SOURCES       := $(filter-out $(SRC_DIR)/main_%.cc,$(ALL_SRC_FILES))
OBJECTS       := $(SOURCES:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

# Main file detection - check both root directory and src directory
# Look for main_coupled
ifneq ($(wildcard main_coupled.cc),)
    MAIN_COUPLED_SRC := main_coupled.cc
    MAIN_COUPLED_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_coupled.cc),)
    MAIN_COUPLED_SRC := main_coupled.cc
    MAIN_COUPLED_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_coupled.cpp),)
    MAIN_COUPLED_SRC := main_coupled.cpp
    MAIN_COUPLED_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_coupled.cpp),)
    MAIN_COUPLED_SRC := main_coupled.cpp
    MAIN_COUPLED_DIR := $(SRC_DIR)
else
    MAIN_COUPLED_SRC :=
    MAIN_COUPLED_DIR :=
endif

# Look for main_russian_doll or main_russian
ifneq ($(wildcard main_russian_doll.cc),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cc
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian_doll.cc),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cc
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian_doll.cpp),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cpp
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian_doll.cpp),)
    MAIN_RUSSIAN_SRC := main_russian_doll.cpp
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian.cc),)
    MAIN_RUSSIAN_SRC := main_russian.cc
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian.cc),)
    MAIN_RUSSIAN_SRC := main_russian.cc
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else ifneq ($(wildcard main_russian.cpp),)
    MAIN_RUSSIAN_SRC := main_russian.cpp
    MAIN_RUSSIAN_DIR := .
else ifneq ($(wildcard $(SRC_DIR)/main_russian.cpp),)
    MAIN_RUSSIAN_SRC := main_russian.cpp
    MAIN_RUSSIAN_DIR := $(SRC_DIR)
else
    MAIN_RUSSIAN_SRC :=
    MAIN_RUSSIAN_DIR :=
endif

# Extract base names for targets
MAIN_COUPLED := $(basename $(MAIN_COUPLED_SRC))
MAIN_RUSSIAN := $(basename $(MAIN_RUSSIAN_SRC))

# Executables in current directory (not in bin/)
TARGET_COUPLED := $(MAIN_COUPLED)
TARGET_RUSSIAN := $(MAIN_RUSSIAN)

# Set default target
ifneq ($(MAIN_COUPLED),)
    DEFAULT_TARGET := coupled
else ifneq ($(MAIN_RUSSIAN),)
    DEFAULT_TARGET := russian
else
    DEFAULT_TARGET := 
endif

# ==============================================================================
# Compiler Flags
# ==============================================================================

# Include paths
INCLUDES      := -I$(INC_DIR) -I$(GETFEM_INC)

# Base compiler flags with system header suppression
CXXFLAGS_BASE := $(CXXSTANDARD) $(INCLUDES) -Wall -Wextra

# Add VERBOSE flag if enabled
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

# Warning suppression flags - comprehensive set
WARN_SUPPRESS := -Wno-comment \
                 -Wno-unused-but-set-variable \
                 -Wno-unused-parameter \
                 -Wno-pragmas \
                 -Wno-deprecated-declarations \
                 -Wno-misleading-indentation

# Combined compiler flags
CXXFLAGS      := $(CXXFLAGS_BASE) $(CXXFLAGS_OPT) $(DEFINES) $(WARN_SUPPRESS)

# Linker flags
LDFLAGS       := -L$(GETFEM_LIB)

# Libraries - proper order matters! 
# GetFEM should come first, then its dependencies
LDLIBS        := -lgetfem

# Add BLAS/LAPACK libraries
# First try to find via pkg-config
LAPACK_LIBS   := $(shell pkg-config --libs lapack 2>/dev/null)
BLAS_LIBS     := $(shell pkg-config --libs blas 2>/dev/null)

# If pkg-config fails, try to detect OpenBLAS or use standard libraries
ifeq ($(LAPACK_LIBS),)
    # Check for OpenBLAS
    ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblas.so*),)
        BLAS_LIBS := -lopenblas
    else ifneq ($(wildcard /usr/lib/x86_64-linux-gnu/libopenblas.so*),)
        BLAS_LIBS := -lopenblas
    else ifneq ($(wildcard /usr/lib/libopenblas.so*),)
        BLAS_LIBS := -lopenblas
    else
        # Fallback to standard BLAS/LAPACK
        BLAS_LIBS := -lblas
    endif
    LAPACK_LIBS := -llapack
endif

# Add SuperLU library (for sparse linear systems)
SUPERLU_LIB := -lsuperlu

# Add Qhull library (for computational geometry)
QHULL_LIB := -lqhull_r

# Combine all libraries in proper order:
# GetFEM -> SuperLU -> LAPACK -> BLAS -> Qhull -> system libs
LDLIBS        += $(SUPERLU_LIB) $(LAPACK_LIBS) $(BLAS_LIBS) $(QHULL_LIB) -lpthread -lm

# ==============================================================================
# Phony Targets
# ==============================================================================

.PHONY: all coupled russian clean distclean help info directories check-getfem list-mains

# Default target
all: $(DEFAULT_TARGET)

# ==============================================================================
# Validation Targets
# ==============================================================================

check-getfem:
	@echo "==> Checking GetFEM installation..."
	@if [ ! -d "$(GETFEM_INC)" ]; then \
		echo "ERROR: GetFEM include directory not found: $(GETFEM_INC)"; \
		echo "Please set GETFEM_PREFIX to your GetFEM installation path:"; \
		echo "  make GETFEM_PREFIX=/path/to/getfem"; \
		exit 1; \
	fi
	@if [ ! -d "$(GETFEM_LIB)" ]; then \
		echo "ERROR: GetFEM library directory not found: $(GETFEM_LIB)"; \
		echo "Please set GETFEM_PREFIX to your GetFEM installation path:"; \
		echo "  make GETFEM_PREFIX=/path/to/getfem"; \
		exit 1; \
	fi
	@echo "==> GetFEM found at: $(GETFEM_PREFIX)"

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
ifneq ($(MAIN_COUPLED),)
	@echo "  make coupled  -> builds $(MAIN_COUPLED)"
endif
ifneq ($(MAIN_RUSSIAN),)
	@echo "  make russian  -> builds $(MAIN_RUSSIAN)"
endif

# ==============================================================================
# Main Build Targets
# ==============================================================================

coupled: check-getfem directories $(TARGET_COUPLED)
	@echo "==> Build complete: $(TARGET_COUPLED)"
	@echo "==> Run with: ./$(TARGET_COUPLED)"
ifeq ($(USE_MPI),1)
	@echo "==> (MPI enabled - use mpirun/mpiexec to run)"
endif

russian: check-getfem directories $(TARGET_RUSSIAN)
	@echo "==> Build complete: $(TARGET_RUSSIAN)"
	@echo "==> Run with: ./$(TARGET_RUSSIAN)"
ifeq ($(USE_MPI),1)
	@echo "==> (MPI enabled - use mpirun/mpiexec to run)"
endif

# ==============================================================================
# Linking Rules
# ==============================================================================

$(TARGET_COUPLED): $(OBJECTS) $(OBJ_DIR)/$(MAIN_COUPLED).o
	@echo "==> Linking $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

$(TARGET_RUSSIAN): $(OBJECTS) $(OBJ_DIR)/$(MAIN_RUSSIAN).o
	@echo "==> Linking $@"
	@$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

# ==============================================================================
# Compilation Rules
# ==============================================================================

# Compile .cc files from src directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 2>&1 | grep -vE "(note: '#pragma message|this is the location)" || true

# Compile .cc files from root directory (for main files)
$(OBJ_DIR)/%.o: %.cc | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 2>&1 | grep -vE "(note: '#pragma message|this is the location)" || true

# Compile .cpp files from src directory
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 2>&1 | grep -vE "(note: '#pragma message|this is the location)" || true

# Compile .cpp files from root directory (for main files)
$(OBJ_DIR)/%.o: %.cpp | $(OBJ_DIR)
	@echo "==> Compiling $<"
	@$(CXX) $(CXXFLAGS) -c $< -o $@ 2>&1 | grep -vE "(note: '#pragma message|this is the location)" || true

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
	@rm -f $(MAIN_COUPLED) $(MAIN_RUSSIAN)

distclean: clean
	@echo "==> Removing all generated files"
	@find . -name "*~" -delete
	@find . -name "*.o" -delete
	@find . -name "core" -delete

info:
	@echo "===================================================================="
	@echo "Build Configuration"
	@echo "===================================================================="
	@echo "Compiler:        $(CXX)"
	@echo "MPI Enabled:     $(USE_MPI)"
	@echo "C++ Standard:    $(CXXSTANDARD)"
	@echo "Build Type:      $(BUILD_TYPE)"
	@echo "Verbose Mode:    $(VERBOSE_MODE)"
	@echo "GETFEM Prefix:   $(GETFEM_PREFIX)"
	@echo "GETFEM Include:  $(GETFEM_INC)"
	@echo "GETFEM Library:  $(GETFEM_LIB)"
	@echo "--------------------------------------------------------------------"
ifneq ($(MAIN_COUPLED_SRC),)
	@echo "Main coupled:    $(MAIN_COUPLED_DIR)/$(MAIN_COUPLED_SRC)"
	@echo "  -> Output:     ./$(MAIN_COUPLED)"
else
	@echo "Main coupled:    NOT FOUND"
endif
ifneq ($(MAIN_RUSSIAN_SRC),)
	@echo "Main russian:    $(MAIN_RUSSIAN_DIR)/$(MAIN_RUSSIAN_SRC)"
	@echo "  -> Output:     ./$(MAIN_RUSSIAN)"
else
	@echo "Main russian:    NOT FOUND"
endif
	@echo "Default target:  $(DEFAULT_TARGET)"
	@echo "--------------------------------------------------------------------"
	@echo "Source files:    $(words $(SOURCES)) files"
	@echo "Object files:    $(words $(OBJECTS)) objects"
	@echo "--------------------------------------------------------------------"
	@echo "Compiler Flags:  $(CXXFLAGS)"
	@echo "Linker Flags:    $(LDFLAGS)"
	@echo "Libraries:       $(LDLIBS)"
	@echo "===================================================================="

help:
	@echo "===================================================================="
	@echo "Available Targets"
	@echo "===================================================================="
	@echo "  make [all]           - Build default target ($(DEFAULT_TARGET))"
ifneq ($(MAIN_COUPLED),)
	@echo "  make coupled         - Build $(MAIN_COUPLED) in current directory"
endif
ifneq ($(MAIN_RUSSIAN),)
	@echo "  make russian         - Build $(MAIN_RUSSIAN) in current directory"
endif
	@echo "  make clean           - Remove object files and binaries"
	@echo "  make distclean       - Remove all generated files"
	@echo "  make info            - Display build configuration"
	@echo "  make help            - Display this help message"
	@echo "  make check-getfem    - Verify GetFEM installation"
	@echo "  make list-mains      - List detected main source files"
	@echo ""
	@echo "===================================================================="
	@echo "Build Options (set via environment or command line)"
	@echo "===================================================================="
	@echo "  BUILD_TYPE=<type>    - Set build type"
	@echo "                         Values: debug, release (default), profile"
	@echo "  VERBOSE_MODE=<0|1>   - Enable VERBOSE preprocessor flag"
	@echo "                         0 = disabled (default), 1 = enabled"
	@echo "  USE_MPI=<0|1>        - Enable MPI parallel compilation"
	@echo "                         0 = g++ (default), 1 = mpic++"
	@echo "  GETFEM_PREFIX=<path> - Set GetFEM installation prefix"
	@echo "                         Default: auto-detected or /usr/local"
	@echo ""
	@echo "===================================================================="
	@echo "Usage Examples"
	@echo "===================================================================="
	@echo "  make                              # Build default target"
	@echo "  make russian                      # Build russian target"
	@echo "  make BUILD_TYPE=debug             # Build with debug symbols"
	@echo "  make VERBOSE_MODE=1               # Build with VERBOSE defined"
	@echo "  make USE_MPI=1                    # Build with MPI support"
	@echo "  make USE_MPI=1 VERBOSE_MODE=1     # MPI + verbose"
	@echo "  make -j4                          # Parallel compilation (4 cores)"
	@echo "  make list-mains                   # See where main files are"
	@echo "  make BUILD_TYPE=debug VERBOSE_MODE=1"
	@echo "  make GETFEM_PREFIX=/custom/path   # Use custom GetFEM location"
	@echo ""
	@echo "  # Without pkg-config:"
	@echo "  make GETFEM_PREFIX=/home/user/getfem/getfem-5.4.4"
	@echo ""
	@echo "===================================================================="
	@echo "Required Libraries"
	@echo "===================================================================="
	@echo "  - GetFEM (finite element library)"
	@echo "  - BLAS/LAPACK or OpenBLAS (linear algebra)"
	@echo "  - SuperLU (sparse linear systems)"
	@echo "  - Qhull (computational geometry)"
	@echo ""
	@echo "  Install on Ubuntu/Debian:"
	@echo "    sudo apt-get install libsuperlu-dev libqhull-dev"
	@echo "    sudo apt-get install libopenblas-dev liblapack-dev"
	@echo "===================================================================="

# ==============================================================================
# Dependency Generation
# ==============================================================================

# Automatically generate dependencies for all object files
DEPS := $(OBJECTS:.o=.d)
ifneq ($(MAIN_COUPLED),)
    DEPS += $(OBJ_DIR)/$(MAIN_COUPLED).d
endif
ifneq ($(MAIN_RUSSIAN),)
    DEPS += $(OBJ_DIR)/$(MAIN_RUSSIAN).d
endif

-include $(DEPS)

# Pattern rules for generating dependency files
$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cc | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: %.cc | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@

$(OBJ_DIR)/%.d: %.cpp | $(OBJ_DIR)
	@$(CXX) $(INCLUDES) -MM -MT '$(OBJ_DIR)/$*.o' $< > $@