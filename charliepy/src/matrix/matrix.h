#ifndef CLP_MATRIXMODULE_H
#define CLP_MATRIXMODULE_H

#ifdef __cplusplus
extern "C" {
#endif

/* This setup is fairly explanatory. Here data points to the underlying data in
 * the matrix, which is an array of objects from FLINT. Strides, as in numpy,
 * allows us to use Python's slice notation to create views of the underyling
 * data. As the matrix may simply be a 'view' of another matrix we have to keep
 * track of the matrix which owns the data; this is controlled by the pointer
 * base. If base is a NULL pointer then this matrix owns its data.
 *
 * We have matrix is a pointer to the
 * underlying FLINT matrix. This is going to be defined over an underlying base
 * ring which is determined by the variable ring. Currently the following rings
 * are supported:
 *
 *      integers = 'Z'
 * 
 * All the heavy lifting for this data type is dealt with by FLINT. */
typedef struct {
    PyObject_HEAD
    char *data;
    char ring;
    Py_ssize_t nrows;
    Py_ssize_t ncols;
    Py_ssize_t rstride;
    Py_ssize_t cstride;
    PyObject *base;
} ClpMatrixObject;

#ifdef __cplusplus
}
#endif
#endif
