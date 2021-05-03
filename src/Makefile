# -----------------------------------------------------------------------------
# Setup...
# -----------------------------------------------------------------------------

# Main executables...
EXC = formod obs2spec brightness climatology jsec2time limb nadir planck time2jsec memoryinfo hash

# Library directories...
LIBDIR =
LIBDIR += -L ../lib/build/lib64
LIBDIR += -L ../lib/build/lib

# Include directories...
INCDIR =
INCDIR += -I ../lib/build/include

# Linking...
STATIC = 1

# Profiling...
#PROF = 1

# -----------------------------------------------------------------------------
# Set flags for GNU compiler...
# -----------------------------------------------------------------------------

# Compilers...
CC = gcc
CU = nvcc

# Check whether a CUDA compiler is available...
ENABLEGPU = $(shell command -v $(CU) 2> /dev/null)

# For reproducability compile the git key (hash) into the code
GITKEY = $(shell git log | head -1 | sed -e 's/commit //g')

# Feature flags...
CFEAFLAGS =
CFEAFLAGS += -fopenmp
CFEAFLAGS += -fno-common
CFEAFLAGS += -fshort-enums
CFEAFLAGS += -fno-omit-frame-pointer
# CFEAFLAGS += -fno-finite-math-only
# CFEAFLAGS += -D DIRECTORY_WITH_GAS_NAME
CFEAFLAGS += -D SHOW_GIT_KEY=$(GITKEY)

# Compiler optimization flags...
COPTFLAGS=
COPTFLAGS += -O3
# COPTFLAGS += -Ofast
# COPTFLAGS += -g
# COPTFLAGS += -pg
# COPTFLAGS += -mtune=power8
# COPTFLAGS += -mcpu=power8
# COPTFLAGS= -O0 -mtune=power8 -mcpu=power8 -g -DGPUDEBUG

CUOPTFLAGS += -D USE_UNIFIED_MEMORY_FOR_TABLES

# Compiler flags...
CFLAGS = -std=gnu99 $(INCDIR) -DHAVE_INLINE -DGSL_DISABLE_DEPRACTED $(COPTFLAGS) $(CERRFLAGS) $(CFEAFLAGS) -pedantic -Winline -Wall -W -Wconversion -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wno-long-long -Wfatal-errors -Wno-deprecated-declarations

# CUDA Compiler optimization flags...
## CUOPTFLAGS = -O3 -gencode arch=compute_60,code=sm_60 $(INCDIR) ## sm_60 for P100
CUOPTFLAGS += -O3 -gencode arch=compute_70,code=sm_70 $(INCDIR)

# CUDA Compiler flags...
CUFLAGS =
CUFLAGS += --std=c++11
CUFLAGS += --use_fast_math $(CUOPTFLAGS) $(patsubst %,-Xcompiler %,$(FLAGS)) 
CUFLAGS += -Xptxas -v  
CUFLAGS += -D RETURN_AFTER_INIT
#CUFLAGS +=  -DGPUDEBUG  

# Linker flags...
LDFLAGS = $(LIBDIR) -lgsl -lgslcblas -lm

# Objects...
OBJECTS = jurassic.o CPUdrivers.o
ifdef ENABLEGPU
	OBJECTS += GPUdrivers.o
	CFEAFLAGS += -DhasGPU
	LDFLAGS += -lcudart
	### module CUDA/10.1.168 in JSC workstation cluster does not provide this path, why?
 	## LIBDIR += -L /usr/local/zam/CUDA/cuda-10.1.168/lib64
 	LIBDIR += -L $(CUDA_PATH)/lib64 
	### furthermore it requires 
	###      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/zam/CUDA/cuda-10.1.168/lib64
	### to run
else
	# Linking...
	ifdef STATIC
		CFLAGS += -static
	endif
endif

# Dependencies...
DEP = jurassic.d CPUdrivers.d GPUdrivers.d $(EXEC:%=%.d)

# -----------------------------------------------------------------------------
# Targets...
# -----------------------------------------------------------------------------

all: $(EXC)

.PHONY: allclean clean
clean:
	rm -f *.o *~ *.d

allclean: clean
	rm -f $(EXC)

bak:
	tar chzf jurassic-gpu_`date +"%y%m%d%H%M"`.tgz Makefile *.c *.h

indent:
	indent -br -brf -brs -bfda -ce -cdw -lp -npcs -npsl *.c *.h

$(EXC): %: %.c $(OBJECTS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Generic rules for all object files...
%.o:%.c %.d Makefile
	$(CC) $(CFLAGS) -c $<
ifdef ENABLEGPU
%.o:%.cu %.d Makefile
	$(CU) $(CUFLAGS) -c $<
else
%.o:%.c %.d Makefile
	$(CC) $(CFLAGS) -c $<
endif

# Dependency rules for correct building in case of changed .h files...
%.d:%.c Makefile
	$(CC) $(CFLAGS) -MM $< -MF $@
%.d:%.cu Makefile
	$(CC) $(CFLAGS) -MM $< -MF $@
-include $(DEP)
