#include <Python.h>
#include <stddef.h>

/* Functions */
static PyObject *
_chartabs_charcolA(PyObject *self, PyObject *args)
{
    PyObject *schm, *t, *col;
    Py_ssize_t k, ind, len_schm;

    if (!PyArg_ParseTuple(args, "OOn", &schm, &t, &k))
        return NULL;

    /* All entries of col are NULL until set by PyList_SetItem */
    len_schm = PyList_GET_SIZE(schm);
    col = PyList_New(len_schm);

    for (ind = 0; ind < len_schm; ind++) {
        PyObject *pi, *val;
        Py_ssize_t j, len_pi;

        val = PyLong_FromLong(0L);

        /* Make pi point to pi[k] */
        pi = PyList_GET_ITEM(schm, ind);
        pi = PyList_GET_ITEM(pi, k);
        len_pi = PyList_GET_SIZE(pi);

        for (j = 0; j < len_pi; j++) {
            PyObject *i, *t_i, *tmp;
            Py_ssize_t i_ssize_t;

            i = PyList_GET_ITEM(pi, j);
            i_ssize_t = PyLong_AsSsize_t(i);

            if (i_ssize_t < 0) {
                i_ssize_t = -i_ssize_t;
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                tmp = PyNumber_Subtract(val, t_i);
            }
            else {
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                tmp = PyNumber_Add(val, t_i);
            }

            Py_DECREF(val);
            val = tmp;
        }
        
        /* Steals a reference to val */ 
        PyList_SET_ITEM(col, ind, val);
    }

    return col;
}

static PyObject *
_chartabs_charcolB(PyObject *self, PyObject *args)
{
    PyObject *schm, *t, *col;
    Py_ssize_t k, len_schm, ind;
    long p, q;

    /* Note p and q are either 0 or 1. */
    if (!PyArg_ParseTuple(args, "OOnll", &schm, &t, &k, &p, &q))
        return NULL;
    
    /* All entries of col are NULL until set by PyList_SetItem */
    len_schm = PyList_GET_SIZE(schm);
    col = PyList_New(len_schm);

    for (ind = 0; ind < len_schm; ind++) {
        PyObject *val, *pi, *pi_k;
        Py_ssize_t j, len_pi_k;

        val = PyLong_FromLong(0L);

        pi = PyList_GET_ITEM(schm, ind);

        /* Make pi_k point to pi[0][k], all references here are borrowed  */
        pi_k = PyList_GET_ITEM(pi, (Py_ssize_t)0);
        pi_k = PyList_GET_ITEM(pi_k, k);
        len_pi_k = PyList_GET_SIZE(pi_k);

        for (j = 0; j < len_pi_k; j++) {
            PyObject *list, *i, *swp, *t_i, *sum;
            long swp_long, i_ssize_t;

            /* Each entry of pi[0][k] is a list [a, b] such that a is
             * the signed index of the value we're interested in and b
             * is 0 or 1. */
            list = PyList_GET_ITEM(pi_k, j);
            i = PyList_GET_ITEM(list, (Py_ssize_t)0);
            swp = PyList_GET_ITEM(list, (Py_ssize_t)1);
            swp_long = PyLong_AsLong(swp);

            i_ssize_t = PyLong_AsSsize_t(i);

            if (i_ssize_t < 0) {
                i_ssize_t = -i_ssize_t;
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                if (swp_long && q)
                    sum = PyNumber_Add(val, t_i);
                else
                    sum = PyNumber_Subtract(val, t_i);
            }
            else {
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                if (swp_long && q)
                    sum = PyNumber_Subtract(val, t_i);
                else 
                    sum = PyNumber_Add(val, t_i);
            }

            Py_DECREF(val);
            val = sum;
        }
        
        /* Make pi_k point to pi[1][k], all references here are borrowed  */
        pi_k = PyList_GET_ITEM(pi, (Py_ssize_t)1);
        pi_k = PyList_GET_ITEM(pi_k, k);
        len_pi_k = PyList_GET_SIZE(pi_k);

        for (j = 0; j < len_pi_k; j++) {
            PyObject *i, *t_i, *sum;
            long i_ssize_t;

            i = PyList_GET_ITEM(pi_k, j);

            i_ssize_t = PyLong_AsSsize_t(i);

            if (i_ssize_t < 0) {
                i_ssize_t = -i_ssize_t;
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                if (p)
                    sum = PyNumber_Add(val, t_i);
                else
                    sum = PyNumber_Subtract(val, t_i);
            }
            else {
                t_i = PyList_GET_ITEM(t, --i_ssize_t);
                if (p)
                    sum = PyNumber_Subtract(val, t_i);
                else 
                    sum = PyNumber_Add(val, t_i);
            }

            Py_DECREF(val);
            val = sum;
        }

        /* Steals a reference to val */ 
        PyList_SET_ITEM(col, ind, val);
    }

    return col;
};

/* Docstrings */
PyDoc_STRVAR(_chartabs__doc__,
"Provides C implementations of functions for computing columns of\n"
"characters tables of type A and B Weyl groups.");

PyDoc_STRVAR(charcolA__doc__,
"Returns a column of the character table of S_{n+1} from a column of\n"
"the character table of S_{k+1}.");

PyDoc_STRVAR(charcolB__doc__,
"Returns a column of the character table of W_n from a column of\n"
"the character table of W_k.");

/* Module method table */
static PyMethodDef
_chartabs_methods[] = {
    {"charcolA", _chartabs_charcolA, METH_VARARGS, charcolA__doc__},
    {"charcolB", _chartabs_charcolB, METH_VARARGS, charcolB__doc__},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/* Module definition */
static struct PyModuleDef
_chartabs = {
    PyModuleDef_HEAD_INIT,
    "_chartabs",
    _chartabs__doc__,
    -1,
    _chartabs_methods
};

/* Initialize the module */
PyMODINIT_FUNC
PyInit__chartabs(void)
{
    PyObject *m = PyModule_Create(&_chartabs);
    if (m == NULL)
        return NULL;
    return m;
}




