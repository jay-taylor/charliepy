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

#define CLP_PERM_MODULE
#include "permutat.h"

/* Constants used for detecting the rank of a permutation. */
static const unsigned long MAX1 = 255;
static const unsigned long MAX2 = 65535;
static const unsigned long MAX4 = 4294967295;

/* We maintain a dynamic buffer on the Python private heap for use
 * throughout this module. This circumvents endless calls to malloc.
 * Any function is free to use the buffer. We note that the buffer is
 * *not* assumed to have any default value. Those functions needing a
 * default value should remember to set it before using the buffer. It
 * is perfectly legitimate to resize the buffer when needed but it is
 * important to remember to update the buffer size.  Here _BUFFSIZE
 * records the size of the memory block pointed to by perm_buff
 * in bytes (where byte means sizeof(char)). We note that the buffer
 * should never be shrunk, this simply just doesn't make sense. */
static size_t _BUFFSIZE;
static char *perm_buff;

/* Utility functions for the buffer. */
static char *
ClpPerm_Buffer(void)
{
    return perm_buff;
}

static int
ClpPerm_ResizeBuffer(size_t newsize)
{
    char *tmp;

    if (newsize > _BUFFSIZE) {
        tmp = PyMem_Resize(perm_buff, char, newsize);
        if (tmp == NULL) {
            PyErr_SetNone(PyExc_MemoryError);
            return -1;
        }
        perm_buff = tmp;
        _BUFFSIZE = newsize;
    }

    return 0;
}

static void
ClpPerm_ZeroBuffer(unsigned long n)
{
    unsigned long i;

    if (n <= MAX1) {
        for (i = 0; i < n; i++) {
            ((ClpPerm_UINT1 *) perm_buff)[i] = 0;
        }
    }
    else if (n <= MAX2) {
        for (i = 0; i < n; i++) {
            ((ClpPerm_UINT2 *) perm_buff)[i] = 0;
        }
    }
    else {
        for (i = 0; i < n; i++) {
            ((ClpPerm_UINT4 *) perm_buff)[i] = 0;
        }
    }
}

/* This macro, taken from the GAP3 source code, returns the image of I
 * under the permutation PT of degree DG. Don't use this with arguments
 * that have side effects! */
#define IMAGE(I, PT, DG) ( (I) < (DG) ? (PT)[(I)] : (I) )

static PyTypeObject ClpPerm_Type; /* Forward */

/* Functions for 8-bit integers. */
#define _PERM_T 1
#include "permutat_templates.c"

/* Functions for 16-bit integers. */
#define _PERM_T 2
#include "permutat_templates.c"

/* Functions for 32-bit integers. */
#define _PERM_T 4
#include "permutat_templates.c"

/* We now must include all the seperate multiplication functions. There
 * are 9 in total. As above the multiplication function is then just a
 * series of switch statements which hands off to the appropriate
 * function. */
#define _PERM_T 1
#define _PERM_S 1
#define _PERM_M 1
#include "permutat_mul_templates.c"

#define _PERM_T 1
#define _PERM_S 2
#define _PERM_M 2
#include "permutat_mul_templates.c"

#define _PERM_T 1
#define _PERM_S 4
#define _PERM_M 4
#include "permutat_mul_templates.c"

#define _PERM_T 2
#define _PERM_S 1
#define _PERM_M 2
#include "permutat_mul_templates.c"

#define _PERM_T 2
#define _PERM_S 2
#define _PERM_M 2
#include "permutat_mul_templates.c"

#define _PERM_T 2
#define _PERM_S 4
#define _PERM_M 4
#include "permutat_mul_templates.c"

#define _PERM_T 4
#define _PERM_S 1
#define _PERM_M 4
#include "permutat_mul_templates.c"

#define _PERM_T 4
#define _PERM_S 2
#define _PERM_M 4
#include "permutat_mul_templates.c"

#define _PERM_T 4
#define _PERM_S 4
#define _PERM_M 4
#include "permutat_mul_templates.c"


/***** INCLUDES OF TEMPLATES MUST OCCUR BEFORE THIS POINT *****/

static PyObject *
perm_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *initial;
    size_t len;

    if (type == &ClpPerm_Type
            && !_PyArg_NoKeywords("permutat.Perm()", kwds)) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O:Perm", &initial)) {
        return NULL;
    }

    /* Check we got a sequence or a permutation. */
    if (ClpPerm_Check(initial)) {
        unsigned short typ = ClpPerm_RANK(initial);

        switch (typ) {
            case 1:
                return permcopy_1(initial);
            case 2:
                return permcopy_2(initial);
            case 4:
                return permcopy_4(initial);
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
ClpPerm_New(unsigned long degree)
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

static PyObject *
ClpPerm_Identity(unsigned long degree)
{
    unsigned long i;
    PyObject *identity;

    if (degree <= MAX1) {
        identity = newperm_NO_CHECKS_1(degree);

        if (identity == NULL) {
            return NULL;
        }

        /* Set the array to the identity permutation. */
        for (i = 0; i < degree; i++) {
            ((ClpPerm_UINT1 *) ClpPerm_ARRAY(identity))[i] = i;
        }
    }
    else if (degree <= MAX2) {
        identity = newperm_NO_CHECKS_2(degree);

        if (identity == NULL) {
            return NULL;
        }

        /* Set the array to the identity permutation. */
        for (i = 0; i < degree; i++) {
            ((ClpPerm_UINT2 *) ClpPerm_ARRAY(identity))[i] = i;
        }
    }
    else {
        identity = newperm_NO_CHECKS_4(degree);

        if (identity == NULL) {
            return NULL;
        }

        /* Set the array to the identity permutation. */
        for (i = 0; i < degree; i++) {
            ((ClpPerm_UINT4 *) ClpPerm_ARRAY(identity))[i] = i;
        }
    }

    return identity;
}


static void
perm_dealloc(ClpPermObject* self)
{
    if (self->array != NULL) {
        PyMem_Del(self->array);
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
perm_repr(ClpPermObject *self)
{
    unsigned short typ = self->rank;

    switch (typ) {
        case 1:
            return perm_repr_1(self);
        case 2:
            return perm_repr_2(self);
        case 4:
            return perm_repr_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_CycleDecomp(ClpPermObject *self)
{
    unsigned short typ = self->rank;

    switch (typ) {
        case 1:
            return perm_cycledecomp_1(self);
        case 2:
            return perm_cycledecomp_2(self);
        case 4:
            return perm_cycledecomp_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_CycleType(ClpPermObject *self, void *closure)
{
    unsigned short typ = self->rank;

    switch (typ) {
        case 1:
            return perm_cycletype_1(self);
        case 2:
            return perm_cycletype_2(self);
        case 4:
            return perm_cycletype_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}


static PyObject *
ClpPerm_CycleStructure(ClpPermObject *self, void *closure)
{
    unsigned short typ = self->rank;

    switch (typ) {
        case 1:
            return perm_cyclestructure_1(self);
        case 2:
            return perm_cyclestructure_2(self);
        case 4:
            return perm_cyclestructure_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_Order(ClpPermObject *self, void *closure)
{
    unsigned short typ = self->rank;

    switch (typ) {
        case 1:
            return perm_order_1(self);
        case 2:
            return perm_order_2(self);
        case 4:
            return perm_order_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_Cycle(PyObject *self, PyObject *args)
{
    PyObject *pt_py;
    unsigned long pt;

    if (!PyArg_ParseTuple(args, "O!:cycle", &PyLong_Type, &pt_py)) {
        return NULL;
    }

    pt = PyLong_AsUnsignedLong(pt_py);
    if (pt >= MAX4 || (pt == (unsigned long) -1 && PyErr_Occurred())) {
        PyErr_Format(PyExc_ValueError,
            "permutations are only defined on integers in the range "
            "[0, ..., %lu]", MAX4 - 1);
        return NULL;
    }

    if (ClpPerm_DEGREE(self) <= MAX1) {
        return cycletuple_1((ClpPermObject *) self, pt);
    }
    else if (ClpPerm_DEGREE(self) <= MAX2) {
        return cycletuple_2((ClpPermObject *) self, pt);
    }
    else {
        return cycletuple_4((ClpPermObject *) self, pt);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_Cycles(PyObject *self, PyObject *args)
{
    if (ClpPerm_DEGREE(self) <= MAX1) {
        return perm_cycles_1((ClpPermObject *) self);
    }
    else if (ClpPerm_DEGREE(self) <= MAX2) {
        return perm_cycles_2((ClpPermObject *) self);
    }
    else {
        return perm_cycles_4((ClpPermObject *) self);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_ToList(PyObject *self, PyObject *args)
{
    if (ClpPerm_DEGREE(self) <= MAX1) {
        return perm_tolist_1((ClpPermObject *) self);
    }
    else if (ClpPerm_DEGREE(self) <= MAX2) {
        return perm_tolist_2((ClpPermObject *) self);
    }
    else {
        return perm_tolist_4((ClpPermObject *) self);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_OnSets(ClpPermObject *self, PyObject *args)
{
    PyObject *list;
    Py_ssize_t len;
    size_t ulen;
    unsigned short typ = ClpPerm_RANK(self);

    if (!PyArg_ParseTuple(args, "O!:onsets", &PyList_Type, &list)) {
        return NULL;
    }

    len = PyList_Size(list);
    if (len == -1) {
        return NULL;
    }

    ulen = (size_t) len;
    if (ulen > MAX4) {
        PyErr_SetString(PyExc_ValueError,
            "The list contains too many sets to define a permutation.");
        return NULL;
    }

    if (ulen <= MAX1) {
        switch (typ) {
            case 1:
                return perm_onsets_11(self, list, ulen);
            case 2:
                return perm_onsets_21(self, list, ulen);
            case 4:
                return perm_onsets_41(self, list, ulen);
        }
    }
    else if (ulen <= MAX2) {
        switch (typ) {
            case 1:
                return perm_onsets_12(self, list, ulen);
            case 2:
                return perm_onsets_22(self, list, ulen);
            case 4:
                return perm_onsets_42(self, list, ulen);
        }
    }
    else {
        switch (typ) {
            case 1:
                return perm_onsets_14(self, list, ulen);
            case 2:
                return perm_onsets_24(self, list, ulen);
            case 4:
                return perm_onsets_44(self, list, ulen);
        }
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_Permute(ClpPermObject *self, PyObject *args)
{
    PyObject *seq;
    unsigned short typ = ClpPerm_RANK(self);

    if (!PyArg_ParseTuple(args, "O:permute", &seq)) {
        return NULL;
    }

    if (seq == NULL) {
        return NULL;
    }

    if (PyList_Check(seq)) {
        switch (typ) {
            case 1:
                return perm_permute_list_1(self, seq);
            case 2:
                return perm_permute_list_2(self, seq);
            case 4:
                return perm_permute_list_4(self, seq);
        }
    }

    if (PyTuple_Check(seq)) {
        switch (typ) {
            case 1:
                return perm_permute_tuple_1(self, seq);
            case 2:
                return perm_permute_tuple_2(self, seq);
            case 4:
                return perm_permute_tuple_4(self, seq);
        }
    }

    /* Can't deal with what we were given. */
    PyErr_SetString(PyExc_TypeError,
        "This function only accepts lists or tuples.");
    return NULL;
}


/* The following hash function is taken from the Python tuple hash
 * function. */
static Py_hash_t
perm_hash(ClpPermObject *self)
{
    unsigned long deg = self->degree;

    if (deg <= MAX1) {
        return perm_hash_1(self);
    }
    else if (deg <= MAX2) {
        return perm_hash_2(self);
    }
    else {
        return perm_hash_4(self);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return -1;
}

/*** Module Functions ***/
/* Note the following function is not very efficient because it cycles
 * through all the elements twice. However this is only likely to be
 * used to make it easier for the user to input permutations, so it's
 * not too important that it is not very performant. */
static PyObject *
ClpPerm_Cyclestoperm(PyObject *self, PyObject *args)
{
    unsigned long max;
    PyObject *max_pylong = PyLong_FromLong(0L);
    Py_ssize_t i, numcycles = PyTuple_Size(args);

    if (max_pylong == NULL) {
        return NULL;
    }

    /* If you didn't get any arguments construct the identity in S_0. */
    if (numcycles == 0) {
        if (PyErr_Occurred()) {
            return NULL;
        }
        else {
            return newperm_NO_CHECKS_1(0);
        }
    }

    /* We start by checking that every element in every cycle that
     * we've been passed is an integer and at the same time find the
     * maximum possible entry. */
    for (i = 0; i < numcycles; i++) {
        Py_ssize_t j, cyclen;
        PyObject *cyc_ob, *cycle = PyTuple_GET_ITEM(args, i);

        if (!PySequence_Check(cycle)) {
            Py_DECREF(max_pylong);
            PyErr_SetString(PyExc_ValueError,
                "This function only accepts sequences.");
            return NULL;
        }

        cyclen = PySequence_Size(cycle);
        for (j = 0; j < cyclen; j++) {
            /* New reference. */
            cyc_ob = PySequence_GetItem(cycle, j);
            if (!PyLong_Check(cyc_ob)) {
                Py_XDECREF(cyc_ob);
                Py_DECREF(max_pylong);
                PyErr_SetString(PyExc_ValueError,
                    "This function only accepts sequences whose entries are integers.");
                return NULL;
            }

            if (PyObject_RichCompareBool(max_pylong, cyc_ob, Py_LT) == 1) {
                Py_DECREF(max_pylong);
                max_pylong = cyc_ob;
            }
            else {
                Py_DECREF(cyc_ob);
            }
        }
    }

    /* Now we can start constructing the permutation. */
    max = PyLong_AsUnsignedLong(max_pylong);
    if (max == (unsigned long) -1 && PyErr_Occurred()) {
        PyErr_Format(PyExc_ValueError,
            "elements of a cycle must lie in the range %lu to %lu",
            0L, MAX4 - 2);
        Py_DECREF(max_pylong);
        return NULL;
    }
    Py_DECREF(max_pylong);

    /* From this point there are no living Python objects. */
    if (max > MAX4 - 2) {
        PyErr_Format(PyExc_ValueError,
            "elements of a cycle must lie in the range %lu to %lu",
            0L, MAX4 - 2);
        return NULL;
    }

    if (max <= MAX1) {
        return cyclestoperm_1(args, max + 1, numcycles);
    }
    else if (max <= MAX2) {
        return cyclestoperm_2(args, max + 1, numcycles);
    }
    else {
        return cyclestoperm_4(args, max + 1, numcycles);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
ClpPerm_RestrictedPerm(PyObject *self, PyObject *args)
{
    PyObject *perm, *list, *max_pylong = PyLong_FromLong(0L);
    Py_ssize_t i, len;
    unsigned long max, deg;
    unsigned short typ;

    if (max_pylong == NULL) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!O:restrictedperm",
                         &ClpPerm_Type, &perm, &list)) {
        return NULL;
    }
    
    if (!PySequence_Check(list)) {
        PyErr_SetString(PyExc_TypeError,
            "The second argument to this function must be a sequence.");
        return NULL;
    }

    len = PySequence_Size(list);
    if (len == -1) {
        return NULL;
    }
    else if (len == 0) {
        return newperm_NO_CHECKS_1(0);
    }

    /* Find the maximum item in the restricted list. */
    for (i = 0; i < len; i++) {
        PyObject *item = PySequence_GetItem(list, i);

        if (item == NULL) {
            return NULL;
        }
        
        if (!PyLong_Check(item)) {
            PyErr_SetString(PyExc_TypeError,
                "The sequence must only contain positive integers.");
            Py_DECREF(item);
            return NULL;
        }

        if (PyObject_RichCompareBool(max_pylong, item, Py_LT) == 1) {
            Py_DECREF(max_pylong);
            max_pylong = item;
        }
        else {
            Py_DECREF(item);
        }
    }

    max = PyLong_AsUnsignedLong(max_pylong);
    if ((max == (unsigned long) -1 && PyErr_Occurred()) || max > MAX4 - 2) {
        Py_DECREF(max_pylong);
        PyErr_Format(PyExc_ValueError,
            "The sequence must only contain integers between %lu and %lu.",
            0L, MAX4 - 2);
        return NULL;
    }
    Py_DECREF(max_pylong);

    deg = max + 1;
    typ = ClpPerm_RANK(perm);

    if (deg <= MAX1) {
        switch (typ) {
            case 1:
                return restrictedperm_11(perm, list, len, deg);
            case 2:
                return restrictedperm_12(perm, list, len, deg);
            case 4:
                return restrictedperm_14(perm, list, len, deg);
        }
    }
    else if (deg <= MAX2) {
        switch (typ) {
            case 1:
                return restrictedperm_21(perm, list, len, deg);
            case 2:
                return restrictedperm_22(perm, list, len, deg);
            case 4:
                return restrictedperm_24(perm, list, len, deg);
        }
    }
    else {
        switch (typ) {
            case 1:
                return restrictedperm_41(perm, list, len, deg);
            case 2:
                return restrictedperm_42(perm, list, len, deg);
            case 4:
                return restrictedperm_44(perm, list, len, deg);
        }
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}




/***** INCLUDES OF MULTIPLICATION TEMPLATES MUST OCCUR BEFORE THIS POINT *****/

/* Note that the variables L, R, typL, typR are not part of the macro
 * definition. It is assumed that they have already been defined. */
#define _BIFUNC_RETURN(X)                                           \
switch (typL) {                                                     \
    case 1: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _11)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _12)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _14)(L, R);     \
        }                                                           \
    }                                                               \
    case 2: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _21)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _22)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _24)(L, R);     \
        }                                                           \
    }                                                               \
    case 4: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _41)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _42)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(perm_, X), _44)(L, R);     \
        }                                                           \
    }                                                               \
}

static PyObject *
perm_mul(PyObject *L, PyObject *R)
{
    unsigned short typL, typR;

    if (ClpPerm_Check(L) && ClpPerm_Check(R)) {
        typL = ClpPerm_RANK(L);
        typR = ClpPerm_RANK(R);

        _BIFUNC_RETURN(mul_perm_perm)
    }
    else if (ClpPerm_Check(L) && PyTuple_Check(R)) {
        PyObject *max_pylong;
        unsigned long max;
        Py_ssize_t i, cyclen = PyTuple_Size(R);

        typL = ClpPerm_RANK(L);
    
        /* If the tuple is empty construct the identity in S_0. */
        if (cyclen == 0 || cyclen == 1) {
            return L;
        }
        else if (cyclen == -1) {
            return NULL;
        }

        max_pylong = PyLong_FromLong(0L);
        if (max_pylong == NULL) {
            return NULL;
        }


        /* We start by checking that every element in the cycle that
         * we've been passed is an integer and at the same time find the
         * maximum possible entry. */
        for (i = 0; i < cyclen; i++) {
            PyObject *cyc_ob = PyTuple_GET_ITEM(R, i);
            if (!PyLong_Check(cyc_ob)) {
                Py_DECREF(max_pylong);
                return NULL;
            }

            if (PyObject_RichCompareBool(max_pylong, cyc_ob, Py_LT) == 1) {
                Py_DECREF(max_pylong);
                Py_INCREF(cyc_ob);
                max_pylong = cyc_ob;
            }
        }

        /* Now we can start constructing the permutation. */
        max = PyLong_AsUnsignedLong(max_pylong);
        if (max == (unsigned long) -1 && PyErr_Occurred()) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX4 - 2);
            Py_DECREF(max_pylong);
            return NULL;
        }
        Py_DECREF(max_pylong);

        /* From this point there are no living Python objects. */
        if (max > MAX4 - 2) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX4 - 2);
            return NULL;
        }

        if (max <= MAX1) {
            R = cycletoperm_1(R, max + 1, cyclen);
            typR = 1;
        }
        else if (max <= MAX2) {
            R = cycletoperm_2(R, max + 1, cyclen);
            typR = 2;
        }
        else {
            R = cycletoperm_4(R, max + 1, cyclen);
            typR = 4;
        }

        if (R == NULL) {
            return NULL;
        }

        _BIFUNC_RETURN(mul_perm_perm)

    }
    else if (PyTuple_Check(L) && ClpPerm_Check(R)) {
        PyObject *max_pylong;
        unsigned long max;
        Py_ssize_t i, cyclen = PyTuple_Size(L);

        typR = ClpPerm_RANK(R);
    
        /* If the tuple is empty construct the identity in S_0. */
        if (cyclen == 0 || cyclen == 1) {
            return L;
        }
        else if (cyclen == -1) {
            return NULL;
        }

        max_pylong = PyLong_FromLong(0L);
        if (max_pylong == NULL) {
            return NULL;
        }

        /* We start by checking that every element in the cycle that
         * we've been passed is an integer and at the same time find the
         * maximum possible entry. */
        for (i = 0; i < cyclen; i++) {
            /* New reference. */
            PyObject *cyc_ob = PyTuple_GET_ITEM(L, i);
            if (!PyLong_Check(cyc_ob)) {
                Py_DECREF(max_pylong);
                return NULL;
            }

            if (PyObject_RichCompareBool(max_pylong, cyc_ob, Py_LT) == 1) {
                Py_DECREF(max_pylong);
                Py_INCREF(cyc_ob);
                max_pylong = cyc_ob;
            }
        }

        /* Now we can start constructing the permutation. */
        max = PyLong_AsUnsignedLong(max_pylong);
        if (max == (unsigned long) -1 && PyErr_Occurred()) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX4 - 2);
            Py_DECREF(max_pylong);
            return NULL;
        }
        Py_DECREF(max_pylong);

        /* From this point there are no living Python objects. */
        if (max > MAX4 - 2) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX4 - 2);
            return NULL;
        }

        if (max <= MAX1) {
            L = cycletoperm_1(L, max + 1, cyclen);
            typL = 1;
        }
        else if (max <= MAX2) {
            L = cycletoperm_2(L, max + 1, cyclen);
            typL = 2;
        }
        else {
            L = cycletoperm_4(L, max + 1, cyclen);
            typL = 4;
        }

        _BIFUNC_RETURN(mul_perm_perm)
    }

    Py_RETURN_NOTIMPLEMENTED;
}

static PyObject *
perm_div(PyObject *L, PyObject *R)
{
    if (ClpPerm_Check(L) && ClpPerm_Check(R)) {
        unsigned short typL = ClpPerm_RANK(L);
        unsigned short typR = ClpPerm_RANK(R);

        _BIFUNC_RETURN(div_perm_perm)
    }
    else if (PyLong_Check(L) && ClpPerm_Check(R)) {
        unsigned short typR = ClpPerm_RANK(R);

        switch (typR) {
            case 1:
                return perm_div_long_perm_1(L, R);
            case 2:
                return perm_div_long_perm_2(L, R);
            case 4:
                return perm_div_long_perm_4(L, R);
        }
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

/* We only define p**i when p is a permutation and i is a long. */
static PyObject *
perm_pow(PyObject *L, PyObject *R, PyObject *Z)
{
    unsigned short typL;

    if (!ClpPerm_Check(L) || !PyLong_Check(R)) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    typL = ClpPerm_RANK(L);

    switch (typL) {
        case 1:
            return perm_pow_1(L, R);
        case 2:
            return perm_pow_2(L, R);
        case 4:
            return perm_pow_4(L, R);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

/* If i is an integer and p, q are permutations then we have
 * 
 *      i^p = the image of i under p,
 *      p^q = q**-1 * p * q,
 *
 * in other words, the latter is simply the conjugate of p by q.*/
static PyObject *
perm_xor(PyObject *L, PyObject *R)
{
    if (ClpPerm_Check(L) && ClpPerm_Check(R)) {
        unsigned short typL = ClpPerm_RANK(L);
        unsigned short typR = ClpPerm_RANK(R);

        _BIFUNC_RETURN(xor_perm_perm)
    }
    else if (PyLong_Check(L) && ClpPerm_Check(R)) {
        unsigned short typR = ClpPerm_RANK(R);

        switch (typR) {
            case 1:
                return perm_xor_long_perm_1(L, R);
            case 2:
                return perm_xor_long_perm_2(L, R);
            case 4:
                return perm_xor_long_perm_4(L, R);
        }
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

/* This function returns True if the permutation is not the identity and
 * False otherwise. This mirrors standards accross Python. */
static int
perm_bool(ClpPermObject *P)
{
    unsigned short typ = P->rank;

    switch (typ) {
        case 1:
            return perm_bool_1(P);
        case 2:
            return perm_bool_2(P);
        case 4:
            return perm_bool_4(P);
    }

    /* If we got here something is wrong. We default to True. */
    PyErr_BadInternalCall();
    return 1;
}

/* Note that the variables c, L, R, typL, typR are not part of the macro
 * definition. It is assumed that they have already been defined. */
#define _BIFUNC_CMP(X)                                              \
switch (typL) {                                                     \
    case 1: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _11)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _12)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _14)(L, R);        \
        }                                                           \
        break;                                                      \
    }                                                               \
    case 2: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _21)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _22)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _24)(L, R);        \
        }                                                           \
        break;                                                      \
    }                                                               \
    case 4: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _41)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _42)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(perm_, X), _44)(L, R);        \
        }                                                           \
    }                                                               \
}

static PyObject *
perm_richcompare(PyObject *L, PyObject *R, int op)
{
    PyObject *res;
    int c;
    unsigned short typL, typR;

    if (!ClpPerm_Check(L) || !ClpPerm_Check(R)) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    typL = ClpPerm_RANK(L);
    typR = ClpPerm_RANK(R);

    switch (op) {
        case Py_LT:
            _BIFUNC_CMP(LT) break;
        case Py_LE:
            _BIFUNC_CMP(LE) break;
        case Py_EQ:
            _BIFUNC_CMP(EQ) break;
        case Py_NE:
            _BIFUNC_CMP(NE) break;
        case Py_GT:
            _BIFUNC_CMP(GT) break;
        case Py_GE:
            _BIFUNC_CMP(GE)
    }

    res = c ? Py_True : Py_False;
    Py_INCREF(res);
    return res;
}

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

/*********************** Permutation Iterator **************************/

static PyTypeObject ClpPermIter_Type; /* Forward */

typedef struct {
    PyObject_HEAD
    unsigned long it_index;
    ClpPermObject *it_perm; /* Set to NULL when iterator is exhausted */
} ClpPermIterObject;

static PyObject *
perm_iter(PyObject *P)
{
    ClpPermIterObject *iter = PyObject_New(ClpPermIterObject,
                                           &ClpPermIter_Type);
    if (iter == NULL) {
        return NULL;
    }
    iter->it_index = 0;
    Py_INCREF(P);
    iter->it_perm = (ClpPermObject *) P;
    return (PyObject *) iter;
}

static void
permiter_dealloc(ClpPermIterObject *iter)
{
    Py_XDECREF(iter->it_perm);
    PyObject_Del(iter);
}

static PyObject *
permiter_next(ClpPermIterObject *iter)
{
    ClpPermObject *P = iter->it_perm;
    unsigned long i = iter->it_index;

    if (P == NULL) {
        return NULL;
    }

    if (i < ClpPerm_DEGREE(P)) {
        PyObject *item;

        switch (ClpPerm_RANK(P)) {
            case 1:
                item = PyLong_FromUnsignedLong(((ClpPerm_UINT1 *) P->array)[i]);
                break;
            case 2:
                item = PyLong_FromUnsignedLong(((ClpPerm_UINT2 *) P->array)[i]);
                break;
            case 4:
                item = PyLong_FromUnsignedLong(((ClpPerm_UINT4 *) P->array)[i]);
        }

        ++iter->it_index;
        return item;
    }

    Py_DECREF(P);
    iter->it_perm = NULL;
    return NULL;
}

static PyTypeObject ClpPermIter_Type = {
    PyVarObject_HEAD_INIT(&PyType_Type, 0)
    "permutat.Perm_Iter",                       /* tp_name */
    sizeof(ClpPermIterObject),                  /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor) permiter_dealloc,              /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_reserved */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                         /* tp_flags */
    0,                                          /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    PyObject_SelfIter,                          /* tp_iter */
    (iternextfunc) permiter_next,               /* tp_iternext */
    0,                                          /* tp_methods */
    0,                                          /* tp_members */
};

/*********************** Module Initialisation **************************/

/* C API methods table */
static _ClpPermAPIMethods _ClpPermAPI = {
    ClpPerm_Buffer,
    ClpPerm_ResizeBuffer,
    ClpPerm_ZeroBuffer,
    ClpPerm_New,
    ClpPerm_Identity,
    ClpPerm_CycleDecomp,
    ClpPerm_Cycle,
    ClpPerm_Order,
    ClpPerm_CycleType,
    ClpPerm_CycleStructure,
    ClpPerm_Cycles,
    ClpPerm_ToList,
    ClpPerm_OnSets,
    ClpPerm_Permute,
    ClpPerm_RestrictedPerm
};

static PyModuleDef permutat = {
    PyModuleDef_HEAD_INIT,
    "permutat",
    "Module giving Python permutation arithmetic.",
    -1,
    permutat_methods, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_permutat(void)
{
    PyObject *m, *api_obj;

    if (PyType_Ready(&ClpPerm_Type) < 0) {
        return NULL;
    }

    m = PyModule_Create(&permutat);
    if (m == NULL) {
        return NULL;
    }

    /* Initialise the buffer. */
    _BUFFSIZE = 1000 * sizeof(ClpPerm_UINT2);
    perm_buff = PyMem_New(char, _BUFFSIZE);
    if (perm_buff == NULL) {
        return PyErr_NoMemory();
    }

    /* Add the permutation object. */
    Py_INCREF(&ClpPerm_Type);
    PyModule_AddObject(m, "Perm", (PyObject *) &ClpPerm_Type);

    /* Add the C API functions. */
    api_obj = PyCapsule_New((void *) &_ClpPermAPI, "permutat._C_API", NULL);
    if (api_obj != NULL) {
        PyModule_AddObject(m, "_C_API", api_obj);
    }

    return m;
}
