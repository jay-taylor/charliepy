/* Ensure that the main file specified enough information. */
#if !defined(_PERM_T) || !defined(_PERM_S) || !defined(_PERM_M)
#  error "Includer has to set _PERM_T and _PERM_S to be either 1, 2 or 4."
#endif

#define _METHOD_NAME(Name) PERM_CAT(PERM_UCAT(Name, _PERM_T), _PERM_S)
#define _UINT_T PERM_CAT(UINT, _PERM_T)
#define _UINT_S PERM_CAT(UINT, _PERM_S) 
#define _UINT_M PERM_CAT(UINT, _PERM_M) 

static PyObject *
_METHOD_NAME(Perm_mul_perm_perm)(PyObject *L, PyObject *R)
{
    PyObject *res;
    _UINT_M *p;
    unsigned long i;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    if (degL <= degR) {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_S)(degR);
        if (res == NULL)
            return NULL;

        p = (_UINT_M *) Perm_array(res);

        for (i = 0; i < degL; i++)
            p[i] = r[l[i]];
        for (i = degL; i < degR; i++)
            p[i] = r[i];
    }
    else {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_T)(degL);
        if (res == NULL)
            return NULL;

        p = (_UINT_M *) Perm_array(res);

        for (i = 0; i < degL; i++)
            p[i] = IMAGE(l[i], r, degR);
    }

    return res;
}

/* Here we define L/R = L*R*(-1). */
static PyObject *
_METHOD_NAME(Perm_div_perm_perm)(PyObject *L, PyObject *R)
{
    PyObject *res;
    _UINT_M *q;
    unsigned long i;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Resize the buffer and invert R into the buffer. */
    if (resize_permbuff((size_t) degR * sizeof(_UINT_S)) != 0) {
        return NULL;
    }
    for (i = 0; i < degR; i++) {
        ((_UINT_S *) permbuff)[r[i]] = i;
    }

    /* Multiply left permutation with the inverse. */
    if (degL <= degR) {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_S)(degR);
        if (res == NULL) {
            return NULL;
        }

        q = (_UINT_M *) Perm_array(res);

        for (i = 0; i < degL; i++) {
            q[i] = ((_UINT_S *) permbuff)[l[i]];
        }
        for (i = degL; i < degR; i++) {
            q[i] = ((_UINT_S *) permbuff)[i];
        }
    }
    else {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_T)(degL);
        if (res == NULL) {
            return NULL;
        }

        q = (_UINT_M *) Perm_array(res);

        for (i = 0; i < degL; i++) {
            q[i] = IMAGE(l[i], ((_UINT_S *) permbuff), degR);
        }
    }

    return res;
}

static PyObject *
_METHOD_NAME(Perm_xor_perm_perm)(PyObject *L, PyObject *R)
{
    PyObject *res;
    _UINT_M *c;
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);
    unsigned long degC = degL < degR ? degR : degL;

    /* Faster if the degrees of the permutations are the same. */
    if (degL == degR) {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_M)(degC);
        if (res == NULL) {
            return NULL;
        }
        c = (_UINT_M *) Perm_array(res);

        for (j = 0; j < degC; j++) {
            c[r[j]] = r[l[j]];
        }
    }

    else {
        res = PERM_UCAT(newperm_NO_CHECKS, _PERM_M)(degC);
        if (res == NULL) {
            return NULL;
        }
        c = (_UINT_M *) Perm_array(res);

        for (j = 0; j < degC; j++) {
            c[IMAGE(j, r, degR)] = IMAGE(IMAGE(j, l, degL), r, degR);
        }
    }

    return res;
}

static int
_METHOD_NAME(Perm_LT)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                if (l[j] < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                if (j < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                if (l[j] < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                if (l[j] < j) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }

    /* They must now be equal. */
    return 0;
}

static int
_METHOD_NAME(Perm_LE)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                if (l[j] < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                if (j < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                if (l[j] < r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                if (l[j] < j) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }

    /* They must now be equal. */
    return 1;
}

static int
_METHOD_NAME(Perm_EQ)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                return 0;
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                return 0;
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                return 0;
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                return 0;
            }
        }
    }

    /* They must now be equal. */
    return 1;
}

static int
_METHOD_NAME(Perm_NE)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                return 1;
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                return 1;
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                return 1;
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                return 1;
            }
        }
    }

    /* They must now be equal. */
    return 0;
}

static int
_METHOD_NAME(Perm_GT)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                if (l[j] > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                if (j > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                if (l[j] > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                if (l[j] > j) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }

    /* They must now be equal. */
    return 0;
}

static int
_METHOD_NAME(Perm_GE)(PyObject *L, PyObject *R)
{
    unsigned long j;
    _UINT_T *l = (_UINT_T *) Perm_array(L);
    _UINT_S *r = (_UINT_S *) Perm_array(R);
    unsigned long degL = Perm_degree(L);
    unsigned long degR = Perm_degree(R);

    /* Find the smallest point where they differ and compare. */
    if (degL <= degR) {
        for (j = 0; j < degL; j++) {
            if (l[j] != r[j]) {
                if (l[j] > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degL; j < degR; j++) {
            if (j != r[j]) {
                if (j > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }
    else {
        for (j = 0; j < degR; j++) {
            if (l[j] != r[j]) {
                if (l[j] > r[j]) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
        for (j = degR; j < degL; j++) {
            if (l[j] != j) {
                if (l[j] > j) {
                    return 1;
                }
                else {
                    return 0;
                }
            }
        }
    }

    /* They must now be equal. */
    return 1;
}

static PyObject *
_METHOD_NAME(restrictedperm)(PyObject *perm, PyObject *list,
                             unsigned long len, unsigned long deg)
{
    PyObject *res;
    unsigned long i;
    
    res = PERM_UCAT(newperm_NO_CHECKS, _PERM_T)(deg);

    if (res == NULL) {
        return NULL;
    }

    /* Start with the identity. */
    for (i = 0; i < deg; i++) {
        ((_UINT_T *) Perm_array(res))[i] = i;
    }

    for (i = 0; i < len; i++) {
        PyObject *item = PySequence_GetItem(list, i);
        unsigned long pt = PyLong_AsUnsignedLong(item);
        Py_DECREF(item);

        ((_UINT_T *) Perm_array(res))[pt] = ((_UINT_S *) Perm_array(perm))[pt];
    }

    /* Finally we check we have a permutation. */
    if (resize_permbuff((size_t) deg * sizeof(_UINT_T)) != 0) {
        Py_DECREF(res);
        return NULL;
    }
    PERM_UCAT(zero_permbuff, _PERM_T)();

    for (i = 0; i < deg; i++) {
        _UINT_T pt = ((_UINT_T *) Perm_array(res))[i];
        if (((_UINT_T *) permbuff)[pt] != 0) {
            PyErr_SetString(PyExc_ValueError,
                "The permutation does not stabilise the sequence.");
            Py_DECREF(res);
            return NULL;
        }
        else {
            ((_UINT_T *) permbuff)[pt] = 1;
        }
    }

    return res;
}

static PyObject *
_METHOD_NAME(Perm_onsets)(Perm *perm, PyObject *list, unsigned long len)
{
    PyObject *res;
    unsigned long i, j;
    _UINT_S *r, *seen = (_UINT_S *) permbuff;

    /* Construct the permutation we will return. */
    res = PERM_UCAT(newperm_NO_CHECKS, _PERM_S)(len);
    if (res == NULL) {
        return NULL;
    }
    r = (_UINT_S *) Perm_array(res);

    /* We use the buffer to keep track of what we have seen. */
    if (resize_permbuff(len * sizeof(_UINT_S)) != 0) {
        Py_DECREF(res);
        return NULL;
    }
    PERM_UCAT(zero_permbuff, _PERM_S)();

    /* First check we got a list of sets. */
    for (i = 0; i < len; i++) {
        PyObject *set = PyList_GET_ITEM(list, i);
        if (!PySet_Check(set)) {
            PyErr_SetString(PyExc_TypeError,
                "Each element of the list must be a set.");
            return NULL;
        }
        if (set == NULL) {
            return NULL;
        }
    }

    /* Now determine the effect of the permutation on the sets. */
    for (i = 0; i < len; i++) {
        PyObject *set = PyList_GET_ITEM(list, i);

        /* Get an element from the set. */
        PyObject *iter = PyObject_GetIter(set);
        if (iter == NULL) {
            Py_DECREF(res);
            return NULL;
        }
        PyObject *first = PyIter_Next(iter);
        Py_DECREF(iter);
        if (first == NULL) {
            Py_DECREF(res);
            return NULL;
        }

        PyObject *image = PERM_UCAT(Perm_xor_long_perm, _PERM_T)(first,
                                                            (PyObject *) perm);
        Py_DECREF(first);
        if (image == NULL) {
            Py_DECREF(res);
            return NULL;
        }

        /* Look for the first set in the list containing the chosen
         * element. */
        for (j = 0; j < len; j++) {
            PyObject *testset = PyList_GET_ITEM(list, j);

            /* If this is it set the permutation. */
            if (PySet_Contains(testset, image)) {
                if (seen[j] != 0) {
                    PyErr_SetString(PyExc_ValueError,
                        "The permutation doesn't permute the list of sets.");
                    Py_DECREF(image);
                    Py_DECREF(res);
                    return NULL;
                }
                r[i] = j;
                break;
            }
            /* If we are at the end and we haven't seen it then the
             * permutation doesn't preserve the list of sets. */
            else if (j == len - 1) {
                PyErr_SetString(PyExc_ValueError,
                    "The permutation doesn't permute the list of sets.");
                Py_DECREF(image);
                Py_DECREF(res);
                return NULL;
            }
        }
        Py_DECREF(image);
    }

    return res;
}










#undef _METHOD_NAME
#undef _UINT_T
#undef _UINT_S
#undef _UINT_M
#undef _PERM_T
#undef _PERM_S
#undef _PERM_M

