#include "temputils.h"

/* Ensure that the main file specified enough information. */
#if !defined(_PERM_T)
#  error "Includer has to set _PERM_T to be either 1, 2 or 4."
#endif

#define _UINT PERM_CAT(ClpPerm_UINT, _PERM_T)
#define _METHOD_NAME(name) PERM_UCAT(name, _PERM_T)
#define _PERMBUFF(I) ((_UINT *) perm_buff)[I]

static PyObject *
_METHOD_NAME(wordtoclassB)(unsigned long n, PyObject *w)
{
    PyObject *perm, *pos, *neg, *pos_tup, *neg_tup, *res;
    Py_ssize_t i, len;
    unsigned long j, N;
    _UINT t, *p;
    _UINT *seen = (_UINT *) ClpPerm_Buffer();

    /* We use the following representation of the Weyl group of type B_n into
     * S_{2n}. If s_0, ..., s_{n-1} are the generators then we have
     * 
     *       s_i -> (i, i+1)(N-i, N-i-1) if i = 0, ..., n-2
     *       s_{n-1} -> (n-1, n)
     * 
     * Here N = 2n-1. This is very similar to the embedding defined in 4.1.18 of
     * [JK81]. To determine the label of the conjugacy class we construct the
     * corresponding element in S_{2n}. We then determine its cycle type,
     * keeping track of whether the cycles are signed or unsigned. */
    N = 2*n - 1;

    perm = ClpPerm_Identity(N + 1);
    if (perm == NULL) {
        return NULL;
    }
    p = (_UINT *) ClpPerm_ARRAY(perm);

    /* We know that we've been passed a list. */
    len = PyList_GET_SIZE(w);

    /* Lists to hold the lengths of signed/unsigned cycles. */
    pos = PyList_New(0);
    neg = PyList_New(0);
    if (pos == NULL || neg == NULL) {
        Py_XDECREF(pos);
        Py_XDECREF(neg);
        Py_DECREF(perm);
        return NULL;
    }

    /* We actually build the reverse word, which is the inverse of
     * our element. However as any element is conjugate to its
     * inverse this won't matter. */
    for (i = len-1; i >= 0; i--) {
        j = PyLong_AsUnsignedLong(PyList_GET_ITEM(w, i));

        /* Check we're inbounds. */
        if (j >= n) {
            PyErr_SetString(PyExc_ValueError,
                "entries in word are out of bounds");
            Py_DECREF(pos);
            Py_DECREF(neg);
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

    /* Now we compute the signed cycle type. We use the buffer to keep track of
     * points that we've already seen. */
    ClpPerm_ZeroBuffer(n);
    
    for (j = 0; j < n; j++) {
        if (!seen[j]) {
            PyObject *tmp;
            unsigned long k, cyc_len = 1;
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
_METHOD_NAME(wordtoclassD)(unsigned long n, PyObject *w)
{
    PyObject *perm, *pos, *neg, *pos_tup, *neg_tup, *res, *pm_str;
    Py_ssize_t i, len;
    unsigned long g, j, k, N;
    unsigned int flag, non_degen;
    int sgn;
    _UINT t, *p;
    _UINT *seen = (_UINT *) ClpPerm_Buffer();

    /* We use the following representation of the Weyl group of type B_n
     * into S_{2n}. If s_0, ..., s_{n-1} are the generators then we have
     * 
     *       s_i -> (i, i+1)(N-i, N-i-1) if i = 0, ..., n-2
     *       s_{n-1} -> (n-2, n)(n-1, n+1)
     * 
     * Here N = 2n-1. This is compatible with our chosen embedding of D_n into
     * B_n and our embedding of B_n into S_{2n}. Our stratergy for deteremining
     * the parameter is the same as the B_n case. We just have to be careful
     * when determining the +/- labelling. */
    N = 2*n - 1;

    perm = ClpPerm_Identity(N + 1);
    if (perm == NULL) {
        return NULL;
    }
    p = (_UINT *) ClpPerm_ARRAY(perm);

    /* We know that we've been passed a list. */
    len = PyList_GET_SIZE(w);

    /* Lists to hold the lengths of signed/unsigned cycles. */
    pos = PyList_New(0);
    neg = PyList_New(0);
    if (pos == NULL || neg == NULL) {
        Py_XDECREF(pos);
        Py_XDECREF(neg);
        Py_DECREF(perm);
        return NULL;
    }

    /* We actually build the reverse word, which is the inverse of
     * our element. However as any element is conjugate to its
     * inverse this won't matter. */
    for (i = len-1; i >= 0; i--) {
        j = PyLong_AsUnsignedLong(PyList_GET_ITEM(w, i));

        /* Check we're inbounds. */
        if (j >= n) {
            PyErr_SetString(PyExc_ValueError,
                "entries in word are out of bounds");
            Py_DECREF(pos);
            Py_DECREF(neg);
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

    /* Now we compute the signed cycle type. We use the buffer to keep track of
     * points that we've already seen. */
    ClpPerm_ZeroBuffer(n);

    /* We need to know if we have a degenerate class. */
    non_degen = 0;
    
    for (j = 0; j < n; j++) {
        if (!seen[j]) {
            PyObject *tmp;
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

            if (neg_test || cyc_len % 2) {
                non_degen = 1;
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

    /* Obtain the bipartition labelling the class in B_n as tuples. */
    if (PyList_Sort(pos) != 0) {
        Py_XDECREF(pos);
        Py_DECREF(neg);
        Py_DECREF(perm);
        return NULL;
    }
    if (PyList_Sort(neg) != 0) {
        Py_DECREF(pos);
        Py_XDECREF(neg);
        Py_DECREF(perm);
        return NULL;
    }

    /* Reversing the lists can't fail. */
    PyList_Reverse(pos);
    PyList_Reverse(neg);

    pos_tup = PyList_AsTuple(pos);
    Py_DECREF(pos);
    if (pos_tup == NULL) {
        Py_DECREF(neg);
        Py_DECREF(perm);
        return NULL;
    }
    neg_tup = PyList_AsTuple(neg);
    Py_DECREF(neg);
    if (neg_tup == NULL) {
        Py_DECREF(pos_tup);
        Py_DECREF(perm);
        return NULL;
    }

    res = PyTuple_New(2);
    if (res == NULL) {
        Py_DECREF(pos_tup);
        Py_DECREF(neg_tup);
        Py_DECREF(perm);
        return NULL;
    }

    if (non_degen) {
        PyTuple_SET_ITEM(res, 0, pos_tup);
        PyTuple_SET_ITEM(res, 1, neg_tup);
        Py_DECREF(perm);
        return res;
    }

    /* We now have a degenerate class and need to decide the +/- label.  For
     * this we start by finding an element x of the Weyl group of type B_n which
     * conjugates p to an element of the + copy of S_n.  Such an element must
     * exist. If x lies in D_n then we have a + element. If x doesn't lie in D_n
     * then s_{n-1}*x does but s_{n-1} sends a + element to a - element, so in
     * that case it's a minus element. As the D_n subgroup is obtained by
     * intersecting the B_n subgroup with the alternating group A_{2n} < S_{2n}
     * we need only look at the parity of the permutation to determine whether
     * it's in D_n or not. */
    flag = 1;
    sgn = 1;
    while (flag) {
        flag = 0;
        for (g = 0; g < n; g++) {
            j = p[g];
            if (j >= n) {
                /* Conjugate p by (j, N-j). First compute (j, N-j)*p. */
                t = p[j];
                p[j] = p[N-j];
                p[N-j] = t;

                /* Now compute p*(j, N-j). */
                for (k = 0; k < N; k++) {
                    if (p[k] == j) {
                        p[k] = N-j;
                    }
                    else if (p[k] == N-j) {
                        p[k] = j;
                    }
                }

                /* Flip the sign of the conjugating element. */
                sgn *= -1;
                flag = 1;
                break;
            }
        }
    }

    if (sgn == 1) {
        pm_str = PyUnicode_FromString("+");
    }
    else {
        pm_str = PyUnicode_FromString("-");
    }

    if (pm_str == NULL) {
        Py_DECREF(perm);
        Py_DECREF(res);
        return NULL;
    }

    PyTuple_SET_ITEM(res, 0, pos_tup);
    Py_DECREF(neg_tup);
    PyTuple_SET_ITEM(res, 1, pm_str);
    Py_DECREF(perm);
    return res;
}

#undef _UINT
#undef _METHOD_NAME
#undef _PERMBUFF
#undef _PERM_T
