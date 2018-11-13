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

# app for PRS solving SOCP
ALP_PD := $(BINDIR)/ALP_PD

CFLAGS := -g -std=c++11 -MMD -w 
#-march=native -O3 -ffast-math -mfpmath=sse -msse4.2


LIB := -L. -lblas -lgfortran -pthread
INC := -I include /usr/local/include/Eigen
DA := ../data/

all: $(ALP_PD)

$(ALP_PD): build/util.o build/test.o build/algebra.o build/primdual.o
	@echo " $(CC) $^ -o $(ALP_PD) $(LIB)"; $(CC) $^ -o $(ALP_PD) $(LIB)
	@echo " $(ALP_PD) is successfully built."
	@printf '%*s' "150" | tr ' ' "-"
	@printf '\n'

# Compile code to objective files
###################################
$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) $(BINDIR)
	@echo " $(CC) $(CFLAGS) $(INC) -c -o $@ $<"; $(CC) $(CFLAGS) $(INC) -c -o $@ $<

run:
	./$(ALP_PD) -nthread 1 -epoch 100000000 -delta 0.5 -alpha 1.0 -block_size_p 50 -block_size_d 25 -check_step 1000 -A $(DA)rdA25x50.mtx -b $(DA)rdb25x50.mtx -c $(DA)rdc25x50.mtx


# clean up the executables and objective fils
##############################################
clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR) $(BINDIR)"; $(RM) -r $(BUILDDIR) $(BINDIR)

-include $(DEPENDENCY)

