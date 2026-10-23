CC=gcc

# DEBUG
#CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -static -g

# RELEASE
CFLAGS=-pedantic -Wall -std=c99 -D_GNU_SOURCE -fPIC -I/opt/homebrew/opt/lapack/include -I/opt/homebrew/opt/openblas/include -I/opt/homebrew/opt/gsl/include -I../../src
LDFLAGS=-L/opt/homebrew/opt/lapack/lib -L/opt/homebrew/opt/openblas/lib -L/opt/homebrew/opt/gsl/lib
LDLIBS=-lm -llapacke -lopenblas -lgsl -lgslcblas
	
traj: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_traj.c ../../src/utils.c ../../src/models.c -L. $(LDLIBS) -o PCR3BP_traj
	@echo ""
	@echo "Problem file compiled successfully."
	./PCR3BP_traj

meas: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_meas.c ../../src/utils.c ../../src/models.c -L. $(LDLIBS) -o PCR3BP_meas
	@echo ""
	@echo "Problem file compiled successfully."
	./PCR3BP_meas

ukf: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_ukf.c ../../src/utils.c ../../src/models.c ../../src/ukf.c -L. $(LDLIBS) -o PCR3BP_ukf
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./PCR3BP_ukf > ./outputs/PCR3BP_ukf_output.log 2>> ./outputs/PCR3BP_ukf_output.log &
#	./PCR3BP_ukf

cukf: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_cukf.c ../../src/utils.c ../../src/models.c ../../src/cukf.c ../../src/pme.c -L. $(LDLIBS) -o PCR3BP_cukf
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./PCR3BP_cukf > ./outputs/PCR3BP_cukf_output.log 2>> ./outputs/PCR3BP_cukf_output.log &
#	./PCR3BP_cukf

cukf_gbees: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_cukf_gbees.c ../../src/utils.c ../../src/models.c ../../src/cukf.c ../../src/pme.c -L. $(LDLIBS) -o PCR3BP_cukf_gbees
	@echo ""
	@echo "Problem file compiled successfully."
	nohup unbuffer ./PCR3BP_cukf_gbees > ./outputs/PCR3BP_hybrid_cukf.log 2>> ./outputs/PCR3BP_hybrid_cukf.log &
#	./PCR3BP_cukf_gbees

post: 
	@echo ""
	@echo "Compiling problem file ..."
	$(CC) ${CFLAGS} ${LDFLAGS} -Wl,-rpath,./ PCR3BP_posterior.c ../../src/utils.c -L. $(LDLIBS) -o PCR3BP_posterior
	@echo ""
	@echo "Problem file compiled successfully."
	./PCR3BP_posterior

clean:
	@echo "Cleaning up local directory ..."
	@-rm -vf PCR3BP_traj
	@-rm -vf PCR3BP_meas
	@-rm -vf PCR3BP_ukf
	@-rm -vf PCR3BP_cukf
	@-rm -vf PCR3BP_cukf_gbees
	@-rm -vf PCR3BP_posterior