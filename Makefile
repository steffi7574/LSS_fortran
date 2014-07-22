FC=gfortran
PETSC_DIR=/home/sguenther/Documents/software/petsc-3.5.0
PETSC_ARCH=${PETSC_DIR}/arch-darwin-c-debug

vanderpol : vanderpol.f 
	$(FC) -I${PETSC_DIR}/include -I${PETSC_ARCH}/include -cpp -o $@ $^


clean :
	rm -f *.o vanderpol


