#include <Python.h>
#include <stddef.h>
#include "permutat.h"

/* Constants used for detecting the rank of a permutation. */
static const unsigned long MAX1 = 255;
static const unsigned long MAX2 = 65535;

/* Functions for 8-bit integers. */
#define _PERM_T 1
#include "_cdata_templates.c"

/* Functions for 16-bit integers. */
#define _PERM_T 2
#include "_cdata_templates.c"

/* Functions for 32-bit integers. */
#define _PERM_T 4
#include "_cdata_templates.c"

/*********************************************************************** 
 *
 *                     Character Table Methods
 *
 ***********************************************************************/
static PyObject *
charcolA(PyObject *self, PyObject *args)
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
charcolB(PyObject *self, PyObject *args)
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

/*********************************************************************** 
 *
 *                   Class Identification Methods
 *
 ***********************************************************************/
static PyObject *
wordtoclassA(PyObject *self, PyObject *args)
{
    PyObject *word;
    Py_ssize_t i, len;
    unsigned long j, n;
    PyObject *perm, *res;

    /* NOTE: This silently overflows if the Python long representing n
     * is too big for an unsigned long. However this shouldn't cause a
     * problem as n should never be very large. */
    if (!PyArg_ParseTuple(args, "kO!", &n, &PyList_Type, &word)) {
        return NULL;
    }
    
    /* We need to work in the symmetric group S_{n+1}. */
    n++;
    perm = ClpPerm_Identity(n);

    if (word == NULL || perm == NULL) {
        return NULL;
    }

    len = PyList_GET_SIZE(word);

    if (n <= MAX1) {
        ClpPerm_UINT1 t, *p = (ClpPerm_UINT1 *) ClpPerm_ARRAY(perm);
    
        /* We traverse the list in reverse to make the multiplication
         * simpler.*/
        for (i = len-1; i >= 0; i--) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                    "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }

            /* Multiply the permutation by (j, j+1). */
            t = p[j];
            p[j] = p[j+1];
            p[j+1] = t;
        }
    }
    else if (n <= MAX2) {
        ClpPerm_UINT2 t, *p = (ClpPerm_UINT2 *) ClpPerm_ARRAY(perm);
    
        /* We traverse the list in reverse to make the multiplication
         * simpler.*/
        for (i = len-1; i >= 0; i--) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }

            /* Multiply the permutation by (j, j+1). */
            t = p[j];
            p[j] = p[j+1];
            p[j+1] = t;
        }
    }
    else {
        ClpPerm_UINT4 t, *p = (ClpPerm_UINT4 *) ClpPerm_ARRAY(perm);
    
        /* We traverse the list in reverse to make the multiplication
         * simpler.*/
        for (i = len-1; i >= 0; i--) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }

            /* Multiply the permutation by (j, j+1). */
            t = p[j];
            p[j] = p[j+1];
            p[j+1] = t;
        }
    }

    res = ClpPerm_CycleType((ClpPermObject *) perm, NULL);
    Py_DECREF(perm);
    return res;
}

static PyObject *
wordtoclassB(PyObject *self, PyObject *args)
{
    PyObject *word;
    unsigned long n;

    /* NOTE: This silently overflows if the Python long representing n
     * is too big for an unsigned long. However this shouldn't cause a
     * problem as n should never be very large. */
    if (!PyArg_ParseTuple(args, "kO!", &n, &PyList_Type, &word)) {
        return NULL;
    }

    if (word == NULL) {
        return NULL;
    }

    if (2*n <= MAX1) {
        return wordtoclassB_1(n, word);
    }
    else if (2*n <= MAX2) {
        return wordtoclassB_2(n, word);
    }
    else {
        return wordtoclassB_4(n, word);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
wordtoclassD(PyObject *self, PyObject *args)
{
    PyObject *word;
    unsigned long n;

    /* NOTE: This silently overflows if the Python long representing n
     * is too big for an unsigned long. However this shouldn't cause a
     * problem as n should never be very large. */
    if (!PyArg_ParseTuple(args, "kO!", &n, &PyList_Type, &word)) {
        return NULL;
    }

    if (word == NULL) {
        return NULL;
    }

    if (2*n <= MAX1) {
        return wordtoclassD_1(n, word);
    }
    else if (2*n <= MAX2) {
        return wordtoclassD_2(n, word);
    }
    else {
        return wordtoclassD_4(n, word);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

#if 0
static PyObject *
wordtoclassB(PyObject *self, PyObject *args)
{
    PyObject *word, *perm, *pos, *neg, *pos_tup, *neg_tup, *res;
    Py_ssize_t i, len;
    unsigned long j, k, n, N;

    /* NOTE: This silently overflows if the Python long representing n
     * is too big for an unsigned long. However this shouldn't cause a
     * problem as n should never be very large. */
    if (!PyArg_ParseTuple(args, "kO!", &n, &PyList_Type, &word)) {
        return NULL;
    }

    if (word == NULL) {
        return NULL;
    }

    /* We use the following representation of the Weyl group of type B_n
     * into S_{2n}. If s_0, ..., s_{n-1} are the generators then we have
     * 
     *       s_i -> (i, i+1)(N-i, N-i-1) if i = 0, ..., n-2
     *       s_{n-1} -> (n-1, n)
     * 
     * Here N = 2n-1. This is very similar to the embedding defined in
     * 4.1.18 of [JK81]. */
    N = 2*n - 1;
    perm = ClpPerm_Identity(N + 1);

    if (perm == NULL) {
        return NULL;
    }

    len = PyList_GET_SIZE(word);

    /* Lists to hold the lengths of signed/unsigned cycles. */
    pos = PyList_New(0);
    neg = PyList_New(0);
    if (pos == NULL || neg == NULL) {
        Py_XDECREF(pos);
        Py_XDECREF(neg);
        Py_DECREF(perm);
    }

    /* We actually build the reverse word, which is the inverse of
     * our element. However as any element is conjugate to its
     * inverse this won't matter. */
    if (2*n <= MAX1) {
        PyObject *tmp;
        ClpPerm_UINT1 t;
        ClpPerm_UINT1 *p = (ClpPerm_UINT1 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT1 *seen = (ClpPerm_UINT1 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);

        /* We traverse the list in reverse to make the multiplication
         * simpler.*/
        for (i = len-1; i >= 0; i--) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                    "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (j, j+1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }
    else if (2*n <= MAX2) {
        PyObject *tmp;
        ClpPerm_UINT2 t;
        ClpPerm_UINT2 *p = (ClpPerm_UINT2 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT2 *seen = (ClpPerm_UINT2 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);

        for (i = 0; i < len; i++) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (j, j+1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }
    else {
        PyObject *tmp;
        ClpPerm_UINT4 t;
        ClpPerm_UINT4 *p = (ClpPerm_UINT4 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT4 *seen = (ClpPerm_UINT4 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);
    
        for (i = 0; i < len; i++) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (j, j+1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }

    /* Sort and return the bipartition as tuples. */
    if (PyList_Sort(pos) != 0) {
        Py_XDECREF(pos);
        Py_DECREF(neg);
        return NULL;
    }
    if (PyList_Sort(neg) != 0) {
        Py_DECREF(pos);
        Py_XDECREF(neg);
        return NULL;
    }

    /* Reversing the lists can't fail. */
    PyList_Reverse(pos);
    PyList_Reverse(neg);

    pos_tup = PyList_AsTuple(pos);
    Py_DECREF(pos);
    if (pos_tup == NULL) {
        Py_DECREF(neg);
        return NULL;
    }
    neg_tup = PyList_AsTuple(neg);
    Py_DECREF(neg);
    if (neg_tup == NULL) {
        Py_DECREF(pos_tup);
        return NULL;
    }

    res = PyTuple_New(2);
    if (res == NULL) {
        return NULL;
    }
    PyTuple_SET_ITEM(res, 0, pos_tup);
    PyTuple_SET_ITEM(res, 1, neg_tup);

    return res;
}

static PyObject *
wordtoclassD(PyObject *self, PyObject *args)
{
    PyObject *word, *perm, *pos, *neg, *pos_tup, *neg_tup, *res;
    Py_ssize_t i, len;
    unsigned long j, k, n, N;

    /* NOTE: This silently overflows if the Python long representing n
     * is too big for an unsigned long. However this shouldn't cause a
     * problem as n should never be very large. */
    if (!PyArg_ParseTuple(args, "kO!", &n, &PyList_Type, &word)) {
        return NULL;
    }

    if (word == NULL) {
        return NULL;
    }

    /* We use the following representation of the Weyl group of type B_n
     * into S_{2n}. If s_0, ..., s_{n-1} are the generators then we have
     * 
     *       s_i -> (i, i+1)(N-i, N-i-1) if i = 0, ..., n-2
     *       s_{n-1} -> (n-2, n)(n-1, n+1)
     * 
     * Here N = 2n-1. This is compatible with our embedding of the type
     * B_n Weyl group. */
    N = 2*n - 1;
    perm = ClpPerm_Identity(N + 1);

    if (perm == NULL) {
        return NULL;
    }

    len = PyList_GET_SIZE(word);

    /* Lists to hold the lengths of signed/unsigned cycles. */
    pos = PyList_New(0);
    neg = PyList_New(0);
    if (pos == NULL || neg == NULL) {
        Py_XDECREF(pos);
        Py_XDECREF(neg);
        Py_DECREF(perm);
    }

    /* We actually build the reverse word, which is the inverse of
     * our element. However as any element is conjugate to its
     * inverse this won't matter. */
    if (N <= MAX1) {
        PyObject *tmp;
        ClpPerm_UINT1 t;
        ClpPerm_UINT1 *p = (ClpPerm_UINT1 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT1 *seen = (ClpPerm_UINT1 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);

        /* We traverse the list in reverse to make the multiplication
         * simpler.*/
        for (i = len-1; i >= 0; i--) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                    "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (n-2, n)(n-1, n+1). */
                t = p[j];
                p[j] = p[j+2];
                p[j+2] = t;

                t = p[j-1];
                p[j-1] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }
    else if (N <= MAX2) {
        PyObject *tmp;
        ClpPerm_UINT2 t;
        ClpPerm_UINT2 *p = (ClpPerm_UINT2 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT2 *seen = (ClpPerm_UINT2 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);
    
        for (i = 0; i < len; i++) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (n-2, n)(n-1, n+1). */
                t = p[j];
                p[j] = p[j+2];
                p[j+2] = t;

                t = p[j-1];
                p[j-1] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }
    else {
        PyObject *tmp;
        ClpPerm_UINT4 t;
        ClpPerm_UINT4 *p = (ClpPerm_UINT4 *) ClpPerm_ARRAY(perm);
        ClpPerm_UINT4 *seen = (ClpPerm_UINT4 *) ClpPerm_Buffer();

        ClpPerm_ZeroBuffer(n);
    
        for (i = 0; i < len; i++) {
            j = PyLong_AsUnsignedLong(PyList_GET_ITEM(word, i));

            /* Check we're inbounds. */
            if (j >= n) {
                PyErr_SetString(PyExc_ValueError,
                        "entries in word are out of bounds");
                Py_DECREF(perm);
                return NULL;
            }
            else if (j == n-1) {
                /* Multiply the permutation by (n-2, n)(n-1, n+1). */
                t = p[j];
                p[j] = p[j+2];
                p[j+2] = t;

                t = p[j-1];
                p[j-1] = p[j+1];
                p[j+1] = t;
            }
            else {
                /* Multiply the permutation by (j, j+1)(N-j, N-j-1). */
                t = p[j];
                p[j] = p[j+1];
                p[j+1] = t;

                t = p[N-j];
                p[N-j] = p[N-j-1];
                p[N-j-1] = t;
            }
        }

        /* Now we compute the signed cycle type. */
        for (j = 0; j < n; j++) {
            if (!seen[j]) {
                unsigned long cyc_len = 1;
                unsigned int neg_test = 0;

                seen[j] = 1;

                /* Get the cycle length and check if it's signed. */
                for (k = p[j]; k != j; k = p[k]) {
                    cyc_len++;
                    if (k == N-j) {
                        neg_test = 1;
                    }
                    if (k < n) {
                        seen[k] = 1;
                    }
                    else {
                        seen[N-k] = 1;
                    }
                }

                /* Add the cycle length to the correct list. */
                if (neg_test) {
                    tmp = PyLong_FromUnsignedLong(cyc_len/2);
                    if (PyList_Append(neg, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
                else {
                    tmp = PyLong_FromUnsignedLong(cyc_len);
                    if (PyList_Append(pos, tmp) != 0) {
                        Py_XDECREF(tmp);
                        Py_DECREF(pos);
                        Py_DECREF(neg);
                        Py_DECREF(perm);
                        return NULL;
                    }
                    Py_DECREF(tmp);
                }
            }
        }

        /* We now no longer need the permutation. */
        Py_DECREF(perm);
    }

    /* Sort and return the bipartition as tuples. */
    if (PyList_Sort(pos) != 0) {
        Py_XDECREF(pos);
        Py_DECREF(neg);
        return NULL;
    }
    if (PyList_Sort(neg) != 0) {
        Py_DECREF(pos);
        Py_XDECREF(neg);
        return NULL;
    }

    /* Reversing the lists can't fail. */
    PyList_Reverse(pos);
    PyList_Reverse(neg);

    pos_tup = PyList_AsTuple(pos);
    Py_DECREF(pos);
    if (pos_tup == NULL) {
        Py_DECREF(neg);
        return NULL;
    }
    neg_tup = PyList_AsTuple(neg);
    Py_DECREF(neg);
    if (neg_tup == NULL) {
        Py_DECREF(pos_tup);
        return NULL;
    }

    res = PyTuple_New(2);
    if (res == NULL) {
        return NULL;
    }
    PyTuple_SET_ITEM(res, 0, pos_tup);
    PyTuple_SET_ITEM(res, 1, neg_tup);

    return res;
}
#endif



/* Docstrings */
PyDoc_STRVAR(_cdata__doc__,
"Provides C implementations of functions for computing columns of\n"
"characters tables of type A and B Weyl groups.");

PyDoc_STRVAR(charcolA__doc__,
"Returns a column of the character table of S_{n+1} from a column of\n"
"the character table of S_{k+1}.");

PyDoc_STRVAR(charcolB__doc__,
"Returns a column of the character table of W_n from a column of\n"
"the character table of W_k.");

PyDoc_STRVAR(wordtoclassA__doc__,
"Returns the class parameter of a word in the standard generators of\n"
"type A Coxeter group. This is simply a tuple giving the cycletype of\n"
"the corresponding permutation in the symmetric group.");

/* Module method table */
static PyMethodDef
_cdata_methods[] = {
    {"charcolA", charcolA, METH_VARARGS, charcolA__doc__},
    {"charcolB", charcolB, METH_VARARGS, charcolB__doc__},
    {"wordtoclassA", wordtoclassA, METH_VARARGS, wordtoclassA__doc__},
    {"wordtoclassB", wordtoclassB, METH_VARARGS, wordtoclassA__doc__},
    {"wordtoclassD", wordtoclassD, METH_VARARGS, wordtoclassA__doc__},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

/* Module definition */
static struct PyModuleDef
_cdata = {
    PyModuleDef_HEAD_INIT,
    "_cdata",
    _cdata__doc__,
    -1,
    _cdata_methods
};

/* Initialize the module */
PyMODINIT_FUNC
PyInit__cdata(void)
{
    PyObject *m = PyModule_Create(&_cdata);
    if (m == NULL)
        return NULL;

    /* Import the permutat module. */
    if (import_permutat() < 0) {
        return NULL;
    }

    return m;
}




