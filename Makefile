SOURCES=network.f95 lstm.f95 softmax.f95 logreg.f95 ftester.f95 

ftester.so: lstm.a
	f2py -c ftester.pyf lstm.a

lstm.a: $(SOURCES)
	gfortran -cpp -std=f2008ts -shared -fPIC -c $^
	ar rvs $@ $(^:%.f95=%.o)

ftester.pyf: ftester.f95
	f2py ftester.f95 -m ftester -h ftester.pyf

a.out: h5create.f95
	gfortran test-hdf5.f95 -I /usr/include /usr/lib/libhdf5_fortran.so.7
