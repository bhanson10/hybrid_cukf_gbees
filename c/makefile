CC=gcc

REB_DIR=/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/rebound-assist/rebound
ifndef REB_DIR # REBOUND is not in default location and REB_DIR is not set
    $(error ASSIST not in the same directory as REBOUND.  To use a custom location, you Must set the REB_DIR environment variable for the path to your rebound directory, e.g., export REB_DIR=/Users/dtamayo/rebound.)
endif
PROBLEMDIR=$(shell basename `dirname \`pwd\``)"/"$(shell basename `pwd`)
include $(REB_DIR)/src/Makefile.defs
ASSIST_DIR=/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/rebound-assist/assist
HOME_DIR=/Users/behanson/Library/CloudStorage/OneDrive-UCSanDiego/UCSD/Research/Utilities/Functions/c

# RELEASE
CFLAGS=-O3 -march=native -fno-math-errno -fassociative-math -fno-signed-zeros -pedantic -Wall -Wno-c23-extensions -std=c99 -D_GNU_SOURCE -fPIC -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/lapack/include -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/gsl/include -I/opt/homebrew/opt/llvm/include -I$(ASSIST_DIR)/src/ -I$(REB_DIR)/src/ -I./src
LDFLAGS=-L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/lapack/lib -L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/gsl/lib -L/opt/homebrew/opt/llvm/lib
LDLIBS=-lm -llapacke -lopenblas -lgsl -lgslcblas -lomp

gd_NEA: librebound.so libassist.so 
	@echo ""
	@echo "Compiling problem file (OpenMP)..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ gauss_detect_NEA.c ./src/utils.c ./src/utils_assist.c ./src/ukf.c ./src/cukf.c -L. $(LDLIBS) -lassist -lrebound -o gauss_detect_NEA
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./gauss_detect_NEA > ./outputs/gauss_detect_NEA_output.log 2>> ./outputs/gauss_detect_NEA_output.log &
#	./gauss_detect_NEA

gd_SatEnc:
	@echo ""
	@echo "Compiling problem file (OpenMP)..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ gauss_detect_SatEnc.c ./src/utils.c ./src/models.c ./src/ukf.c ./src/cukf.c -L. $(LDLIBS) -o gauss_detect_SatEnc
	@echo ""
	@echo "Problem file compiled successfully."
#	nohup unbuffer ./gauss_detect_SatEnc > ./outputs/gauss_detect_SatEnc_output.log 2>> ./outputs/gauss_detect_SatEnc_output.log &
	./gauss_detect_SatEnc

librebound.so:
	@echo "Compiling shared library librebound.so ..."
	$(MAKE) -C $(REB_DIR)/src/
	@echo "Creating link for shared library librebound.so ..."
	@-rm -f librebound.so
	@ln -s $(REB_DIR)/src/librebound.so .

libassist.so: librebound.so $(ASSIST_DIR)/src/*.c  $(ASSIST_DIR)/src/*.h
	@echo "Compiling shared library libassist.so ..."
	$(MAKE) -C $(ASSIST_DIR)/src/
	@-rm -f libassist.so
	@ln -s $(ASSIST_DIR)/src/libassist.so .

clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf gauss_detect_NEA
	@-rm -vf gauss_detect_SatEnc
	@-rm -vf librebound.so
	@-rm -vf libassist.so