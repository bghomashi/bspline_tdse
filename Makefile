IDIR =  inc/
SDIR =  src/
ODIR =	obj/
DDIR =	dep/
csrc =  $(wildcard src/*.c)

ccsrc = $(wildcard src/*.cpp) \
		$(wildcard src/core/maths/*.cpp) \
		$(wildcard src/core/bspline/*.cpp) \
		$(wildcard src/core/tdse/*.cpp) \
		$(wildcard src/core/tise/*.cpp) \
		$(wildcard src/core/utility/*.cpp) \
		$(wildcard src/eigen_solvers/*.cpp) \
		$(wildcard src/math_libs/*.cpp) \
		$(wildcard src/math_libs/petsc/*.cpp) \
		$(wildcard src/tdse_propagators/*.cpp) \
		$(wildcard src/observables/*.cpp) \
		$(wildcard src/potentials/*.cpp) \
		$(wildcard src/input_validation/*.cpp) \
		$(wildcard src/tests/*.cpp)

OBJ = 	$(subst $(SDIR), $(ODIR), $(csrc:.c=.o)) \
		$(subst $(SDIR), $(ODIR), $(ccsrc:.cpp=.o))
DEP =	$(subst $(ODIR), $(DDIR), $(OBJ:.o=.d))


#preprocessor flags
EXTRAFLAGS = -g -O0 
#CPPFLAGS =  -DGLEW_STATIC -Iinc -IGL -I../freetype2 -I/usr/local/include/freetype2
CPPFLAGS =      -DHAVE_LAPACK_CONFIG_H -DLAPACK_COMPLEX_CPP \
 				-Dlapack_complex_float="std::complex<float>" \
 				-Dlapack_complex_double="std::complex<double>" \
				-D__STDCPP_WANT_MATH_SPEC_FUNCS__ \
				-I$(PETSC_DIR)/include -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
				-I$(SLEPC_DIR)/include -I$(SLEPC_DIR)/$(PETSC_ARCH)/include \
				-Iinc -I/usr/include \
				-Isrc/core -Isrc/core/bspline -Isrc/core/maths -Isrc/core/tdse -Isrc/core/tise \
				-Isrc -std=c++14
#linker flags (directories)
#LDFLAGS
#linker libraries
LDLIBS =-lm -pthread -lhdf5 -L$(SLEPC_DIR)/$(PETSC_ARCH)/lib -L$(PETSC_DIR)/$(PETSC_ARCH)/lib -lpetsc -lslepc
#LDLIBS =-lm -lGL -lGLU -lglfw -lX11 -ldl -lfreetype -pthread -llapacke -lblas

CC = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx
CPP = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx
CXX = $(PETSC_DIR)/$(PETSC_ARCH)/bin/mpicxx

bspline_tdse.out: $(OBJ)
	$(CXX) $(EXTRAFLAGS) -o $@ $^ $(LDLIBS)

-include $(DEP)							#include all dep files in the makefile

.PHONY: clean
clean:
	rm -f $(OBJ) plotter_test
	rm -f $(DEP)


%/ : 
	@mkdir -p $@

.SECONDEXPANSION:
$(DDIR)%.d: $(SDIR)%.cpp | $$(dir $$@)
	@$(CPP) $(CPPFLAGS) $(EXTRAFLAGS) $< -MM -MT $(subst $(DDIR), $(ODIR), $(@:.d=.o)) >$@


$(ODIR)%.o : $(SDIR)%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) $(EXTRAFLAGS) -c $< -o $@

$(ODIR)%.o : ../%.c | $$(dir $$@)
	$(CC) $(CPPFLAGS) -c $< -o $@

$(ODIR)%.o : $(SDIR)%.cpp | $$(dir $$@)
	$(CXX) $(CPPFLAGS) $(EXTRAFLAGS) -c $< -o $@


