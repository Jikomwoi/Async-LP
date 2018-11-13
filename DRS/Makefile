# Makefile for TMAC_SOCP

CC := g++ # This is the main compiler
# CC := clang --analyze # and comment out the linker last line for sanity
# directory of objective files
BUILDDIR := build
# directory of binary files
BINDIR := bin
# directory of source code
SRCDIR := src
# extension of source file
SRCEXT := cc
SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

DEPENDENCY := $(shell find $(BUILDDIR) -type f -name *.d 2>/dev/null)

# app for DRS solving LP
ALP_DRS := $(BINDIR)/ALP_DRS

CFLAGS := -g -std=c++11 -MMD -w 
#-march=native -O3 -ffast-math -mfpmath=sse -msse4.2


LIB := -lblas -lgfortran -pthread 
#LIB := -pthread -lblas -lgfortran -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lgomp 
INC := -I include -I /usr/local/include/Eigen
#INC := -I include -I include/Eigen -I ${MKLROOT}/include
DA := ../data/

all: $(ALP_DRS)

$(ALP_DRS): build/util.o build/test.o build/algebra.o build/DRS.o
	@echo " $(CC) $^ -o $(ALP_DRS) $(LIB)"; $(CC) $^ -o $(ALP_DRS) $(LIB)
	@echo " $(ALP_DRS) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

# Compile code to objective files
###################################
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) $(BINDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -O2 -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -O2 -o $@ $<

run:
	./$(ALP_DRS) -nthread 2 -epoch 1000000 -block_size 25 -check_step 100 -A $(DA)rdA25x50.mtx -b $(DA)rdb25x50.mtx -c $(DA)rdc25x50.mtx -eta 0.5 -lambda 1.0 -scale 1 
# clean up the executables and objective fils
##############################################
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

-include $(DEPENDENCY)

