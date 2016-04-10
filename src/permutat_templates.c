#include "temputils.h"

/* Ensure that the main file specified enough information. */
#if !defined(_PERM_T)
#  error "Includer has to set _PERM_T to be either 1, 2 or 4."
#endif

#define _UINT PERM_CAT(UINT, _PERM_T)
#define _METHOD_NAME(Name) PERM_UCAT(Name, _PERM_T)
#define _PERMBUFF(I) ((_UINT *) permbuff)[I]

/* This is a utility function for the buffer which sets all entries to 0
 * assuming the underling data type of the buffer is _UINT.*/
static void
_METHOD_NAME(zero_permbuff)()
{
    size_t i, len = _BUFFSIZE / sizeof(_UINT);

    for (i = 0; i < len; i++) {
        _PERMBUFF(i) = 0;
    }
}

static int
_METHOD_NAME(setpermitem)(Perm *perm, unsigned long i, PyObject *x, _UINT *seen)
{
    /* Here seen is an array such that seen[i] is 0 if we haven't
     * already come accross i in the list and 1 if we have. */
    unsigned long k;

    if (!PyLong_Check(x)) {
        PyErr_Format(PyExc_ValueError,
            "list element must be an integer in the range %lu to %lu",
            0L, perm->degree - 1);
        return -1;
    }

    k = PyLong_AsUnsignedLong(x);
    /* Check for overflow. */
    if (k == (unsigned long) -1 && PyErr_Occurred()) {
        PyErr_Format(PyExc_ValueError,
            "list element must lie in range %lu to %lu", 0L, perm->degree - 1);
        return -1;
    }

    /* No overflow so can safely check k is in the right range. */
    if (k >= perm->degree || k > MAX_4) {
        PyErr_Format(PyExc_ValueError,
            "list element must lie in range %lu to %lu", 0L, perm->degree - 1);
        return -1;
    }
    
    /* Check if we've already seen it. */
    if (seen[k] != 0) {
        PyErr_Format(PyExc_ValueError,
            "%lu must only occur once in the list", k);
        return -1;
    }
    seen[k] = 1;

    ((_UINT *) perm->array)[i] = k;
    return 0;
}

/* Create a new instance of Perm whose array is set by the list pointed
 * to by initial. There is a seperate function which initiates a new
 * permutation without initiating the underyling array. */
static PyObject *
_METHOD_NAME(newperm)(PyTypeObject *type, unsigned long degree,
                      PyObject *initial)
{
    unsigned long i;
    Perm *self = (Perm *) type->tp_alloc(type, 0);

    if (self == NULL) {
        return NULL;
    }

    if (degree == 0) {
        self->array = NULL;
        self->type = _PERM_T;
        self->degree = 0;
        return (PyObject *) self;
    }

    self->array = PyMem_New(char, (size_t) degree * sizeof(_UINT));
    if (self->array == NULL) {
        Py_DECREF(self);
        return PyErr_NoMemory();
    }
    self->type = _PERM_T;
    self->degree = degree;

    /* Now set the array to the initial list. We need to use the buffer
     * here to keep track of which elements we've already seen, this is
     * passed to setpermitem. */
    if (resize_permbuff((size_t) degree * sizeof(_UINT)) != 0) {
        Py_DECREF(self);
        return NULL;
    }
    _METHOD_NAME(zero_permbuff)();

    /* Set the array elements. */
    for (i = 0; i < degree; i++) {
        /* x is a new reference. */
        PyObject *x = PySequence_GetItem(initial, i);
        if (x == NULL) {
            Py_DECREF(self);
            return NULL;
        }

        if (_METHOD_NAME(setpermitem)(self, i, x, (_UINT *) permbuff) != 0) {
            Py_DECREF(self);
            Py_DECREF(x);
            return NULL;
        }

        Py_DECREF(x);
    }

    return (PyObject *) self;
}

/* This is an internal version of the above function which does not error
 * check the degree and does not initialise the underyling array of the
 * permutation. */
static PyObject *
_METHOD_NAME(newperm_NO_CHECKS)(unsigned long degree)
{
    Perm *self = (Perm *) Perm_Type.tp_alloc(&Perm_Type, 0);

    if (self == NULL) {
        return NULL;
    }

    if (degree == 0) {
        self->array = NULL;
        self->type = _PERM_T;
        self->degree = degree;
    }
    else {
        self->array = PyMem_New(char, (size_t) degree * sizeof(_UINT));
        if (self->array == NULL) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->type = _PERM_T;
        self->degree = degree;
    }

    return (PyObject *) self;
}

/* This is an internal version of the above function which does not error
 * check the degree and does not initialise the underyling array of the
 * permutation. */
static PyObject *
_METHOD_NAME(permcopy)(PyObject *initial)
{
    unsigned long i, deg = Perm_degree(initial);
    PyObject *res = _METHOD_NAME(newperm_NO_CHECKS)(deg);
    _UINT *r, *p;

    if (res == NULL) {
        return NULL;
    }

    r = (_UINT *) Perm_array(res);
    p = (_UINT *) Perm_array(initial);

    for (i = 0; i < deg; i++) {
        r[i] = p[i];
    }

    return res;
}


static PyObject *
_METHOD_NAME(cycletuple)(Perm *self, unsigned long pt)
{
    PyObject *cycle;
    _UINT *p = (_UINT *) self->array;
    _UINT i, j, min, len, count;

    if (pt >= self->degree) {
        /* Get the tuple and fill it. */
        PyObject *cycle = PyTuple_New((Py_ssize_t) 1);
        if (cycle == NULL) {
            return NULL;
        }
    
        PyTuple_SET_ITEM(cycle, 0, PyLong_FromSize_t((size_t) pt));
        return cycle;
    }

    /* Get the length of the cycle and the minimal element. */
    len = 1;
    min = pt;
    for (i = p[pt]; i != pt; i = p[i]) {
        len++;
        if (i < min) {
            min = i;
        }
    }

    /* Get the tuple and fill it. */
    cycle = PyTuple_New((Py_ssize_t) len);
    if (cycle == NULL) {
        return NULL;
    }
    PyTuple_SET_ITEM(cycle, 0, PyLong_FromSize_t((size_t) min));

    count = 1;
    for (j = p[min]; j != min; j = p[j]) {
        PyTuple_SET_ITEM(cycle, count, PyLong_FromSize_t((size_t) j));
        count++;
    }

    return cycle;
}


/* There should be a better way to write the following function that
 * doesn't involve constructing Python objects as intermediaries.
 * However as it's mostly window dressing we'll leave it as is.  */
static PyObject *
_METHOD_NAME(Perm_repr)(Perm *self)
{
    PyObject *s, *blank, *cycles;
    unsigned long deg = self->degree;
    _UINT *p = (_UINT *) self->array;
    _UINT i;

    /* Special case for degree 0. */
    if (deg == 0) {
        return PyUnicode_FromString("()");
    }

    cycles = PyList_New(0L);
    if (cycles == NULL) {
        return NULL;
    }

    for (i = 0; i < deg; i++) {
        if (i != p[i]) {
            _UINT j = p[i];

            /* Only get a cycle if i is the smallest element in the cycle. */
            while (i < j)
                j = p[j];
            if (i == j) {
                PyObject *cycle = _METHOD_NAME(cycletuple)(self, i);
                if (cycle == NULL) {
                    Py_DECREF(cycles);
                    return NULL;
                }
                PyObject *scycle = PyUnicode_FromFormat("%R", cycle);
                if (scycle == NULL) {
                    Py_DECREF(cycles);
                    return NULL;
                }

                Py_DECREF(cycle);
                PyList_Append(cycles, scycle);
                Py_DECREF(scycle);
            }
        }
    }

    /* Special case for the identity. */
    if (PyList_GET_SIZE(cycles) == 0) {
        Py_DECREF(cycles);
        return PyUnicode_FromString("()");
    }
    
    blank = PyUnicode_FromString("");
    s = PyUnicode_Join(blank, cycles);
    Py_DECREF(cycles);
    Py_DECREF(blank);
    return s;
}

static PyObject *
_METHOD_NAME(Perm_cycledecomp)(Perm *self)
{
    PyObject *cycles;
    _UINT *p = (_UINT *) self->array;
    unsigned long i, deg = self->degree;
    
    cycles = PyList_New(0L);
    if (cycles == NULL) {
        return NULL;
    }

    for (i = 0; i < deg; i++) {
        if (p[i] != i) {
            _UINT j = p[i];

            /* Only get a cycle if i is the smallest element in the cycle. */
            while (i < j)
                j = p[j];
            if (i == j) {
                PyObject *cycle = _METHOD_NAME(cycletuple)(self, i);
                if (cycle == NULL) {
                    Py_DECREF(cycles);
                    return NULL;
                }
                if (PyList_Append(cycles, cycle) == -1) {
                    Py_DECREF(cycles);
                    Py_DECREF(cycle);
                    return NULL;
                }
                Py_DECREF(cycle);
            }
        }
    }

    return cycles;
}

static PyObject *
_METHOD_NAME(Perm_cycles)(Perm *self)
{
    PyObject *cycles;
    _UINT *p = (_UINT *) self->array;
    unsigned long i, deg = self->degree;
    
    cycles = PyList_New(0L);
    if (cycles == NULL) {
        return NULL;
    }

    for (i = 0; i < deg; i++) {
        if (p[i] != i) {
            _UINT j = p[i];

            /* Only get a cycle if i is the smallest element in the cycle. */
            while (i < j)
                j = p[j];
            if (i == j) {
                PyObject *cycle = _METHOD_NAME(cycletuple)(self, i);
                if (cycle == NULL) {
                    Py_DECREF(cycles);
                    return NULL;
                }
                if (PyList_Append(cycles, cycle) == -1) {
                    Py_DECREF(cycles);
                    Py_DECREF(cycle);
                    return NULL;
                }
                Py_DECREF(cycle);
            }
        }
        else {
            /* Get the tuple and fill it. */
            PyObject *cycle = PyTuple_New((Py_ssize_t) 1);
            if (cycle == NULL) {
                return NULL;
            }
        
            PyTuple_SET_ITEM(cycle, 0, PyLong_FromUnsignedLong(i));
            if (PyList_Append(cycles, cycle) == -1) {
                Py_DECREF(cycles);
                Py_DECREF(cycle);
                return NULL;
            }
            Py_DECREF(cycle);
        }
    }

    return cycles;
}

static PyObject *
_METHOD_NAME(Perm_cycletype)(Perm *self)
{
    PyObject *res;
    _UINT *cycletype, *p = (_UINT *) self->array;
    unsigned long i, k, len, deg = self->degree;
    
    cycletype = PyMem_New(_UINT, deg + 1);
    if (cycletype == NULL) {
        return NULL;
    }

    /* Initialise the list to 0. */
    for (i = 0; i < deg + 1; i++) {
        cycletype[i] = 0;
    }

    for (i = 0; i < deg; i++) {
        if (p[i] != i) {
            _UINT j = p[i];

            /* Only get a cycle if i is the smallest element in the cycle. */
            while (i < j)
                j = p[j];
            if (i == j) {
                /* Get the length of the cycle. */
                len = 1;
                for (k = p[i]; k != i; k = p[k]) {
                    len++;
                }
                cycletype[len]++;
            }
        }
        else {
            cycletype[1]++;
        }
    }

    res = PyTuple_New(deg + 1);
    if (res == NULL) {
        PyMem_Del(cycletype);
        return NULL;
    }

    for (i = 0; i < deg + 1; i++) {
        PyTuple_SET_ITEM(res, i, PyLong_FromUnsignedLong(cycletype[i]));
    }

    PyMem_Del(cycletype);

    return res;
}

static PyObject *
_METHOD_NAME(Perm_cycletype_part)(Perm *self)
{
    PyObject *res;
    _UINT *cycletype, *p = (_UINT *) self->array;
    unsigned long i, k, len, count, cycnum = 0, deg = self->degree;
    
    cycletype = PyMem_New(_UINT, deg + 1);
    if (cycletype == NULL) {
        return NULL;
    }

    /* Initialise the list to 0. */
    for (i = 0; i < deg + 1; i++) {
        cycletype[i] = 0;
    }

    for (i = 0; i < deg; i++) {
        if (p[i] != i) {
            _UINT j = p[i];

            /* Only get a cycle if i is the smallest element in the cycle. */
            while (i < j)
                j = p[j];
            if (i == j) {
                /* Get the length of the cycle. */
                len = 1;
                for (k = p[i]; k != i; k = p[k]) {
                    len++;
                }
                cycletype[len]++;
                cycnum++;
            }
        }
        else {
            cycletype[1]++;
            cycnum++;
        }
    }

    res = PyTuple_New(cycnum);
    if (res == NULL) {
        PyMem_Del(cycletype);
        return NULL;
    }

    count = 0;

    /* We don't need to check the first entry in the array because it is
     * always equal to 0. */
    for (i = deg; i > 0; i--) {
        for (k = 0; k < cycletype[i]; k++) { 
            PyTuple_SET_ITEM(res, count, PyLong_FromUnsignedLong(i));
            count++;
        }
    }

    PyMem_Del(cycletype);

    return res;
}

static PyObject *
_METHOD_NAME(Perm_tolist)(Perm *self)
{
    unsigned long i, deg = self->degree;
    _UINT *p = (_UINT *) self->array;
    PyObject *res = PyList_New(deg);

    if (res == NULL) {
        return NULL;
    }

    for (i = 0; i < deg; i++) {
        PyList_SET_ITEM(res, i, PyLong_FromUnsignedLong(p[i]));
    }

    return res;
}

static Py_hash_t
_METHOD_NAME(permhash)(Perm *self)
{
    Py_uhash_t x;  /* Unsigned for defined overflow behavior. */
    Py_uhash_t mult = _PyHASH_MULTIPLIER;
    _UINT *p = (_UINT *) self->array;
    unsigned long i, deg = self->degree;

    x = 0x345678UL;
    for (i = 0; i < deg; i++) {
        x = (x ^ p[i]) * mult;
        mult += (Py_hash_t)(82520UL + deg + deg);
    }
    x += 97531UL;
    if (x == (Py_uhash_t)-1)
        x = -2;
    return x;
}


/* The following algorithm works by looping over all cycles and
 * computing the gcd of their lengths. We use the buffer to record when
 * we've seen a value. Note that we first have to make sure the buffer
 * is big enough and is initialised to 0. */
static PyObject *
_METHOD_NAME(Perm_order)(Perm *self)
{
    _UINT i, ord = 1;
    _UINT *p = (_UINT *) self->array;
    unsigned long deg = self->degree;

    if (deg == 0)
        return PyLong_FromLong(1L);

    if (resize_permbuff((size_t) deg * sizeof(_UINT)) != 0)
        return NULL;
    _METHOD_NAME(zero_permbuff)();

    for (i = 0; i < deg; i++) {
        if (_PERMBUFF(i) == 0 && p[i] != i) {
            _UINT j, gcd, s, t, len = 1;

            for (j = p[i]; j != i; j = p[j]) {
                len++;
                _PERMBUFF(j) = 1;
            }

            gcd = len;
            s = ord % len;
            while (s != 0) {
                t = s;
                s = gcd % s;
                gcd = t;
            }
            ord = ord * (len / gcd);
        }
    }

    return PyLong_FromSize_t((size_t) ord);
}

static PyObject *
_METHOD_NAME(cyclestoperm)(PyObject *cycles, unsigned long deg,
                           Py_ssize_t numcycles)
{
    unsigned long i;
    Py_ssize_t j;
    _UINT *seen = (_UINT *) permbuff;
    PyObject *res = _METHOD_NAME(newperm_NO_CHECKS)(deg);
    if (res == NULL) {
        return NULL;
    }

    /* Start by setting the array to the identity permutation. */
    for (i = 0; i < deg; i++) {
        ((_UINT *) Perm_array(res))[i] = i;
    }

    /* We use the buffer to keep track of what we've seen. */
    if (resize_permbuff((size_t) deg * sizeof(_UINT)) != 0) {
        return NULL;
    }
    _METHOD_NAME(zero_permbuff)();

    for (j = 0; j < numcycles; j++) {
        Py_ssize_t k;
        unsigned long m, n;
        PyObject *x, *y;
        PyObject *cycle = PyTuple_GET_ITEM(cycles, j);
        Py_ssize_t cyclen = PySequence_Size(cycle);

        if (cyclen == 0) {
            continue;
        }

        /* Note we know that every element in the cycle is less than deg
         * so they are all in bounds. GetItem gives a new reference. */
        x = PySequence_GetItem(cycle, 0);
        m = PyLong_AsUnsignedLong(x);
        Py_DECREF(x);

        y = PySequence_GetItem(cycle, cyclen-1);
        n = PyLong_AsUnsignedLong(y);
        Py_DECREF(y);

        if (seen[m] != 0) {
            PyErr_Format(PyExc_ValueError,
                "%lu must only occur in one cycle", m);
            Py_DECREF(res);
            return NULL;
        }
        else {
            ((_UINT *) Perm_array(res))[n] = m;
            seen[m] = 1;
        }

        for (k = 1; k < cyclen; k++) {
            y = PySequence_GetItem(cycle, k);
            n = PyLong_AsUnsignedLong(y);
            Py_DECREF(y);

            if (seen[n] != 0) {
                PyErr_Format(PyExc_ValueError,
                    "%lu must only occur in one cycle", n);
                Py_DECREF(res);
                return NULL;
            }
            else {
                ((_UINT *) Perm_array(res))[m] = n;
                seen[n] = 1;
            }
            
            m = n;
        }
    }

    return res;
}

static PyObject *
_METHOD_NAME(cycletoperm)(PyObject *cycle, unsigned long deg,
                          Py_ssize_t cyclen)
{
    unsigned long i, m, n;
    Py_ssize_t k;
    _UINT *seen = (_UINT *) permbuff;
    PyObject *res = _METHOD_NAME(newperm_NO_CHECKS)(deg);
    if (res == NULL) {
        return NULL;
    }

    /* Start by setting the array to the identity permutation. */
    for (i = 0; i < deg; i++) {
        ((_UINT *) Perm_array(res))[i] = i;
    }

    /* We use the buffer to keep track of what we've seen. */
    if (resize_permbuff((size_t) deg * sizeof(_UINT)) != 0) {
        return NULL;
    }
    _METHOD_NAME(zero_permbuff)();

    /* Note we know that every element in the cycle is less than deg
     * so they are all in bounds. */
    m = PyLong_AsUnsignedLong(PyTuple_GET_ITEM(cycle, 0));
    n = PyLong_AsUnsignedLong(PyTuple_GET_ITEM(cycle, cyclen-1));

    ((_UINT *) Perm_array(res))[n] = m;
    seen[m] = 1;

    for (k = 1; k < cyclen; k++) {
        n = PyLong_AsUnsignedLong(PyTuple_GET_ITEM(cycle, k));

        if (seen[n] != 0) {
            PyErr_Format(PyExc_ValueError,
                "%lu must only occur once in the cycle", n);
            Py_DECREF(res);
            return NULL;
        }
        else {
            ((_UINT *) Perm_array(res))[m] = n;
            seen[n] = 1;
        }
        
        m = n;
    }

    return res;
}



/***** NUMERIC METHODS *****/
/* Here we have i/p is the preimage of i under p. */
static PyObject *
_METHOD_NAME(Perm_div_long_perm)(PyObject *L, PyObject *R)
{
    unsigned long i, pre;
    _UINT *r = (_UINT *) Perm_array(R);
    unsigned long degR = Perm_degree(R);

    i = PyLong_AsUnsignedLong(L);
    /* Check for overflow. */
    if (i == (unsigned long) -1 && PyErr_Occurred()) {
        PyErr_Format(PyExc_ValueError,
            "permutation is only defined on integers in the "
            "range [%lu, ..., %lu]", 0L, MAX_4 - 1);
        return NULL;
    }
    
    /* Compute the cycle containing i then the last point hit before
     * coming back to i is the preimage. */
    pre = i;
    if (i < degR) {
        while (r[pre] != i)
            pre = r[pre];
    }

    return PyLong_FromUnsignedLong(pre);
}

/* We only define p**i when p is a permutation and i is a long. */
static PyObject *
_METHOD_NAME(Perm_pow)(PyObject *L, PyObject *R)
{
    PyObject *res;
    long i;
    _UINT *p, *l = (_UINT *) Perm_array(L);
    unsigned long degL = Perm_degree(L);

    i = PyLong_AsLong(R);
    /* Check for overflow. */
    if (i == -1 && PyErr_Occurred()) {
        /* LONG_MAX and LONG_MIN from limits.h included by Python.h. */
        PyErr_Format(PyExc_ValueError,
            "cannot take power of permutation larger than %li or smaller than",
            LONG_MAX, LONG_MIN);
        return NULL;
    }

    res = _METHOD_NAME(newperm_NO_CHECKS)(degL);
    if (res == NULL) {
        return NULL;
    }
    p = (_UINT *) Perm_array(res);

    /* For small powers just compute the power by repeated mapping. */
    if (0 <= i && i < 8) {
        _UINT j, k;
        long e;

        for (j = 0; j < degL; j++) {
            k = j;
            for (e = 0; e < i; e++) {
                k = l[k];
            }
            p[j] = k;
        }
    }

    /* For larger powers we compute the power of the individual cycles. */
    else if (i >= 8) {
        _UINT j, k, s;
        long e, len;

        /* We need to use the buffer to record when we've seen a value. */
        if (resize_permbuff((size_t) degL * sizeof(_UINT)) != 0) {
            Py_DECREF(res);
            return NULL;
        }
        _METHOD_NAME(zero_permbuff)();

        /* Loop over the cycles. */
        for (j = 0; j < degL; j++) {
            if (_PERMBUFF(j) == 0) {
                /* Get the cycle length. */
                len = 1;
                for (k = l[j]; k != j; k = l[k]) {
                    len++;
                    _PERMBUFF(k) = 1;
                }

                /* Raise the cycle to the power i modulo the length. */
                s = j;
                for (e = 0; e < i % len; e++) {
                    s = l[s];
                }
                p[j] = s;

                s = l[s];
                for (k = l[j]; k != j; k = l[k]) {
                    p[k] = s;
                    s = l[s];
                }
            }
        }
    }

    /* The inverse case is easy. */
    else if (i == -1) {
        _UINT j;

        for (j = 0; j < degL; j++) {
            p[l[j]] = j;
        }
    }

    /* As for small positive powers just compute via repeated mapping. */
    else if (i > -8) {
        _UINT j, k;
        long e;

        for (j = 0; j < degL; j++) {
            k = j;
            for (e = 0; e > i; e--) {
                k = l[k];
            }
            p[k] = j;
        }
    }

    /* As for large positive powers we compute the power via cycles. */
    else {
        _UINT j, k, s;
        long e, len;

        /* We need to use the buffer to record when we've seen a value. */
        if (resize_permbuff((size_t) degL * sizeof(_UINT)) != 0) {
            Py_DECREF(res);
            return NULL;
        }
        _METHOD_NAME(zero_permbuff)();

        /* Loop over the cycles. */
        for (j = 0; j < degL; j++) {
            if (_PERMBUFF(j) == 0) {
                /* Get the cycle length. */
                len = 1;
                for (k = l[j]; k != j; k = l[k]) {
                    len++;
                    _PERMBUFF(k) = 1;
                }

                /* Raise the cycle to the power i modulo the length. */
                s = j;
                for (e = 0; e > i % len; e--) {
                    s = l[s];
                }
                p[s] = j;

                s = l[s];
                for (k = l[j]; k != j; k = l[k]) {
                    p[s] = k;
                    s = l[s];
                }
            }
        }
    }

    return res;
}

static PyObject *
_METHOD_NAME(Perm_xor_long_perm)(PyObject *L, PyObject *R)
{
    PyObject *res;
    unsigned long i;

    i = PyLong_AsUnsignedLong(L);
    /* Check for overflow. */
    if (i == (unsigned long) -1 && PyErr_Occurred()) {
        /* LONG_MAX from limits.h included by Python.h. */
        PyErr_Format(PyExc_ValueError,
            "permutations are only defined on integers in the range "
            "[0, ..., %lu]", MAX_4 - 1);
        return NULL;
    }

    if (i >= MAX_4) {
        PyErr_Format(PyExc_ValueError,
            "permutations are only defined on integers in the range "
            "[0, ..., %lu]",
            MAX_4 - 1);
        return NULL;
    }

    /* If we get a NULL pointer no need to DECREF, we'll just pass
     * the error through. */
    if (i >= Perm_degree(R)) {
        res = PyLong_FromLong(i);
    }
    else {
        _UINT tmp = ((_UINT *) Perm_array(R))[i];
        res = PyLong_FromUnsignedLong((unsigned long) tmp);
    }

    return res;
}



#undef _UINT
#undef _METHOD_NAME
#undef _PERMBUFF
#undef _PERM_T

