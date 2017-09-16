/***********************************************************************
 * This file is adapted from the corresponding file, also named
 * permutat.c, in the GAP3 source code. This source code may be obtained
 * from the GAP site
 *
 *      http://www.gap-system.org/Gap3/gap3.html
 *
 * The algorithms are almost always unchanged, except for minor changes
 * in names of variables. However, there are several changes in the
 * implementation vs the GAP implementation. Here we make the
 * distinction between permutations on less than 2**8-1, 2**16-1, and
 * 2**32-1 points. In GAP only the latter 2 distinctions are made. Due
 * to the many different versions of functions that must necessarily be
 * written we have therefore taken a template approach and used heavily
 * the C preprocessor.
 *
 * Of course, the source also looks significantly different due to the
 * inclusion of Python C-API code. The Python specific parts of this
 * module are based on the array module found in the Python standard
 * library.
************************************************************************/
#include <Python.h>
#include <structmember.h>
#include <fmpz.h>

#define CLP_MATRIX_MODULE
#include "matrix.h"

static PyTypeObject ClpMatrix_Type; /* Forward */

static fmpz *
fmpz_mat_entry(ClpMatrixObject *mat, Py_ssize_t i, Py_ssize_t j)
{
    return ((fmpz *) mat->data) + (i*mat->rstride + j*mat->cstride);
}

static fmpz *
pylong_to_fmpz(PyObject *pylong)
{
}


/***** INCLUDES OF TEMPLATES MUST OCCUR BEFORE THIS POINT *****/

static PyObject *
matrix_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *initial;
    size_t len;

    if (type == &ClpMatrix_Type
            && !_PyArg_NoKeywords("matrix.Matrix()", kwds)) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O:Matrix", &initial)) {
        return NULL;
    }

    /* Check we got a sequence or a permutation. */
    if (ClpMatrix_Check(initial)) {
        unsigned short typ = ClpPerm_RANK(initial);

        switch (typ) {
            case 1:
                return matrixcopy_1(initial);
            case 2:
                return matrixcopy_2(initial);
            case 4:
                return matrixcopy_4(initial);
        }
    }
    else if (PySequence_Check(initial)) {
        len = (size_t) PySequence_Size(initial);
    }
    else {
        PyErr_SetString(PyExc_ValueError,
            "permutation must be set by a permutation or sequence");
        return NULL;
    }

    if (len > MAX4) {
        PyErr_SetString(PyExc_ValueError,
            "permutation exceeded maximum size");
        return NULL;
    }

    /* Now get the appropriate type of permutation. */
    if (len <= MAX1) {
        return newperm_1(type, len, initial);
    }
    else if (len <= MAX2) {
        return newperm_2(type, len, initial);
    }
    else {
        return newperm_4(type, len, initial);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpMatrix_New(unsigned long degree)
{
    if (degree <= MAX1) {
        return newperm_NO_CHECKS_1(degree);
    }
    else if (degree <= MAX2) {
        return newperm_NO_CHECKS_2(degree);
    }
    else {
        return newperm_NO_CHECKS_4(degree);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}


static void
matrix_dealloc(ClpPermObject* self)
{
    if (self->array != NULL) {
        PyMem_Del(self->array);
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}

/*** Module Functions ***/

static PyNumberMethods perm_as_number = {
     0,                            /* nb_add */
     0,                            /* nb_subtract */
     (binaryfunc) perm_mul,        /* nb_multiply */
     0,                            /* nb_remainder */
     0,                            /* nb_divmod */
     (ternaryfunc) perm_pow,       /* nb_power */
     0,                            /* nb_negative */
     0,                            /* nb_positive */
     0,                            /* nb_absolute */
     (inquiry) perm_bool,          /* nb_bool */
     0,                            /* nb_invert */
     0,                            /* nb_lshift */
     0,                            /* nb_rshift */
     0,                            /* nb_and */
     (binaryfunc) perm_xor,        /* nb_xor */
     0,                            /* nb_or */
     0,                            /* nb_int */
     0,                            /* *nb_reserved */
     0,                            /* nb_float */
     0,                            /* nb_inplace_add */
     0,                            /* nb_inplace_subtract */
     0,                            /* nb_inplace_multiply */
     0,                            /* nb_inplace_remainder */
     0,                            /* nb_inplace_power */
     0,                            /* nb_inplace_lshift */
     0,                            /* nb_inplace_rshift */
     0,                            /* nb_inplace_and */
     0,                            /* nb_inplace_xor */
     0,                            /* nb_inplace_or */
     (binaryfunc) perm_div,        /* nb_floor_divide */
     (binaryfunc) perm_div,        /* nb_true_divide */
     0,                            /* nb_inplace_floor_divide */
     0,                            /* nb_inplace_true_divide */
     0,                            /* nb_index */
     0,                            /* nb_matrix_multiply */
     0                             /* nb_inplace_matrix_multiply */
};


static PyMethodDef perm_methods[] = {
    {"cycledecomp", (PyCFunction)ClpPerm_CycleDecomp, METH_NOARGS,
     "Give the decomposition of the permutation into disjoint cycles."},
    {"cycle", (PyCFunction)ClpPerm_Cycle, METH_VARARGS,
     "Return the cycle of the permutation which moves the given point."},
    {"cycles", (PyCFunction)ClpPerm_Cycles, METH_NOARGS,
     "Return a list containing all the cycles of the permutation."},
    {"tolist", (PyCFunction)ClpPerm_ToList, METH_NOARGS,
     "Return the permutation as a Python list."},
    {"onsets", (PyCFunction)ClpPerm_OnSets, METH_VARARGS,
     "Return the induced permutation acting on a list of sets."},
    {"permute", (PyCFunction)ClpPerm_Permute, METH_VARARGS,
     "Return the sequence obtained by permuting its entries."},
    {NULL} /* Sentinel */
};

static PyMemberDef perm_members[] = {
    {"degree", T_ULONG, offsetof(ClpPermObject, degree), READONLY,
     "The minimum n such that S_n contains the permutation."},
    {NULL} /* Sentinel */
};

static PyGetSetDef perm_getsetters[] = {
    {"order", (getter)ClpPerm_Order, NULL,
     "The multiplicative order of the permutation."},
    {"cycletype", (getter)ClpPerm_CycleType, NULL,
     "The cycle type of the permutation as a partition."},
    {"cyclestructure", (getter)ClpPerm_CycleStructure, NULL,
     "The lengths of non-trivial cycles of the permutation."},
    {NULL} /* Sentinel */
};


static PyMethodDef permutat_methods[] = {
    {"cyclestoperm", (PyCFunction)ClpPerm_Cyclestoperm, METH_VARARGS,
     "Return a permutation which is a product of the corresponding "
     "disjoint cycles."},
    {"restrictedperm", (PyCFunction)ClpPerm_RestrictedPerm, METH_VARARGS,
     "Return the restriction of a permutation to a subset."},
    {NULL} /* Sentinel */
};

static PyObject *perm_iter(PyObject *); /* Forward */

static PyTypeObject ClpPerm_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "permutat.Perm",                /* tp_name */
    sizeof(ClpPermObject),          /* tp_basicsize */
    0,                              /* tp_itemsize */
    (destructor) perm_dealloc,      /* tp_dealloc */
    0,                              /* tp_print */
    0,                              /* tp_getattr */
    0,                              /* tp_setattr */
    0,                              /* tp_reserved */
    (reprfunc) perm_repr,           /* tp_repr */
    &perm_as_number,                /* tp_as_number */
    0,                              /* tp_as_sequence */
    0,                              /* tp_as_mapping */
    (hashfunc) perm_hash,           /* tp_hash  */
    0,                              /* tp_call */
    0,                              /* tp_str */
    0,                              /* tp_getattro */
    0,                              /* tp_setattro */
    0,                              /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,        /* tp_flags */
    "ClpPerm objects",              /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    perm_richcompare,               /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    perm_iter,                      /* tp_iter */
    0,                              /* tp_iternext */
    perm_methods,                   /* tp_methods */
    perm_members,                   /* tp_members */
    perm_getsetters,                /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    0,                              /* tp_init */
    0,                              /* tp_alloc */
    perm_new,                       /* tp_new */
};

/*********************** Module Initialisation **************************/

static PyModuleDef matrix = {
    PyModuleDef_HEAD_INIT,
    "matrix",
    "Module giving Python access to FLINT matrices.",
    -1,
    matrix_methods, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_matrix(void)
{
    PyObject *m, *api_obj;

    if (PyType_Ready(&ClpMatrix_Type) < 0) {
        return NULL;
    }

    m = PyModule_Create(&matrix);
    if (m == NULL) {
        return NULL;
    }

    /* Add the permutation object. */
    Py_INCREF(&ClpMatrix_Type);
    PyModule_AddObject(m, "Matrix", (PyObject *) &ClpMatrix_Type);

    return m;
}
