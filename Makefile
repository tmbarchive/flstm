a.out: h5create.f95
	gfortran test-hdf5.f95 -I /usr/include /usr/lib/libhdf5_fortran.so.7
