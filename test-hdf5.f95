!
! The following example shows how to create an empty dataset. 
! It creates a file called 'dsetf.h5', defines the
! dataset dataspace, creates a dataset which is a 4x6 integer array,
! and then closes the dataspace, the dataset, and the file.
!

PROGRAM DSETEXAMPLE
USE HDF5 ! This module contains all necessary modules 
IMPLICIT NONE
CHARACTER(LEN=8), PARAMETER :: filename = "dsetf.h5" ! File name
CHARACTER(LEN=4), PARAMETER :: dsetname = "dset"     ! Dataset name
INTEGER(HID_T) :: file_id       ! File identifier 
INTEGER(HID_T) :: dset_id       ! Dataset identifier 
INTEGER(HID_T) :: dspace_id     ! Dataspace identifier
INTEGER(HSIZE_T), DIMENSION(2) :: dims = (/4,6/) ! Dataset dimensions
INTEGER     ::   rank = 2                        ! Dataset rank
INTEGER     ::   error ! Error flag
INTEGER, DIMENSION(4,6) :: dset_data, data_out ! Data buffers
INTEGER(HSIZE_T), DIMENSION(2) :: data_dims
INTEGER i,j
DO i = 1, 4
 DO j = 1, 6
    dset_data(i,j) = (i-1)*6 + j
 END DO
END DO
CALL h5open_f(error)
CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
CALL h5screate_simple_f(rank, dims, dspace_id, error)
CALL h5dcreate_f(file_id, dsetname, H5T_NATIVE_INTEGER, dspace_id, dset_id, error)

data_dims(1) = 4
data_dims(2) = 6
CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dset_data, data_dims, error)

CALL h5dclose_f(dset_id, error)
CALL h5sclose_f(dspace_id, error)
CALL h5fclose_f(file_id, error)
CALL h5close_f(error)
END PROGRAM DSETEXAMPLE 
