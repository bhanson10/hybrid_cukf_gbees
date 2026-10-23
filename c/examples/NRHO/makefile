CC=gcc

# RELEASE
CFLAGS=-O3 -march=native -pedantic -Wall -Wno-c23-extensions -std=c99 -D_GNU_SOURCE -fPIC -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -I/opt/homebrew/opt/lapack/include -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/gsl/include -I../../src
LDFLAGS=-L/opt/homebrew/opt/libomp/lib -L/opt/homebrew/opt/lapack/lib -L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/gsl/lib
LDLIBS=-lm -llapacke -lopenblas -lomp -lgsl -lgslcblas

traj: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ nrho_traj.c ../../src/utils.c ../../src/models.c -L. $(LDLIBS) -o nrho_traj
	@echo ""
	@echo "Problem file compiled successfully."
	./nrho_traj

pf: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ nrho_pf.c ../../src/utils.c ../../src/models.c ../../src/pf.c  -L. $(LDLIBS) -o nrho_pf
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./nrho_pf > ./outputs/nrho_pf_output.log 2>> ./outputs/nrho_pf_output.log &
#	./nrho_pf

cukf: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ nrho_cukf.c ../../src/utils.c ../../src/models.c ../../src/cukf.c -L. $(LDLIBS) -o nrho_cukf
	@echo ""
	@echo "Problem file compiled successfully."
	./nrho_cukf

cukfpf: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ nrho_cukfpf.c ../../src/utils.c ../../src/models.c ../../src/cukf.c ../../src/pf.c -L. $(LDLIBS) -o nrho_cukfpf
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./nrho_cukfpf > ./outputs/nrho_cukfpf_output.log 2>> ./outputs/nrho_cukfpf_output.log &
#	./capstone_ukfpf

clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf nrho_traj
	@-rm -vf nrho_pf
	@-rm -vf nrho_cukf
	@-rm -vf nrho_cukfpf