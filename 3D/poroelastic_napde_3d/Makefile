# ===============================================================================
# A Makefile for compiling a C++ project with sources in 'src/' and headers
# in 'include/', with a main file 'main_coupled.cc' in the project root.
# ===============================================================================

# Compiler
CXX = g++
# CXX = mpic++

# Optimization flags
OPTFLAGS = -O3 -g

# GetFEM path (change to your installation path)
GETFEM_PATH = /home/pacs/getfem/getfem-5.4.4

# Include directories
INCLUDE = -I$(GETFEM_PATH)/include -I$(GETFEM_PATH)/src -I$(GETFEM_PATH)/src/gmm -I./include -I/usr/include

# Preprocessor definitions
DEF_TAGS = -DHAVE_CONFIG -DGMM_USES_BLAS

# Compiler flags (only used for compilation)
CXXFLAGS = -Wall -Wextra $(INCLUDE) $(OPTFLAGS) -MMD -MP $(DEF_TAGS)

# Linker flags and libraries (only used for linking)
LDFLAGS = -rdynamic \
          $(GETFEM_PATH)/src/.libs/libgetfem.a \
          -lsuperlu -llapack -lblas -ldl \
          /usr/lib/x86_64-linux-gnu/libqhull.so.8.0 \
          /usr/lib/x86_64-linux-gnu/liblapack.so.3 \
          /usr/lib/x86_64-linux-gnu/libblas.so.3 \
          /usr/lib/x86_64-linux-gnu/libsuperlu.so

# Executable name
EXEC = main

# Sources folder
FOLDER = src/

# Find all C++ source files in src and add the main file
SRC_FILES = $(wildcard $(FOLDER)*.cc) ./main_coupled.cc
#SRC_FILES = $(wildcard $(FOLDER)*.cc) ./main_russian_doll.cc

# Object files
OBJS = $(SRC_FILES:.cc=.o)

# ===============================================================================
# Rules
# ===============================================================================

.PHONY: all clean

# Main build rule
all: $(EXEC)

# Linking rule
$(EXEC): $(OBJS)
	@echo "Linking executable..."
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
	@echo "Build successful: $(EXEC)"

# Compilation rule
%.o: %.cc
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean rule
clean:
	@echo "Cleaning up build files..."
	@rm -f $(EXEC) $(OBJS) $(OBJS:.o=.d)
	@echo "Clean complete."

# Include dependency files
-include $(OBJS:.o=.d)
