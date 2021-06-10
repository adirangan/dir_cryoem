LIBPATH=/opt/intel/compilers_and_libraries/mac/mkl/lib
LIBS= \
  $(LIBPATH)/libmkl_lapack95_lp64.a \
  $(LIBPATH)/libmkl_blas95_lp64.a \
  $(LIBPATH)/libmkl_intel_lp64.a \
  $(LIBPATH)/libmkl_core.a  \
  $(LIBPATH)/libmkl_sequential.a \
  -lpthread


sptest.ex: sparse_test.f
	ifort -o sptest.ex sparse_test.f  $(LIBS)


