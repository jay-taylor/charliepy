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
#include <stdint.h>
#include "temputils.h"

/* Our permutations will be defined on at most 2**32-2 points. This
 * means that the degree of the permutation is at most 2**32-1, which is
 * the maximum value for an unsigned 32-bit integer. To save memory we
 * consider three different types of permutations, depending upon
 * whether the permutation acts on 2**8-1, 2**16-1, or 2**32-1 points.
 * For this we introduce the following types. */
typedef uint_least8_t UINT1;
typedef uint_least16_t UINT2;
typedef uint_least32_t UINT4;
static const unsigned long MAX_1 = 255;
static const unsigned long MAX_2 = 65535;
static const unsigned long MAX_4 = 4294967295;

/* We maintain a dynamic buffer on the Python private heap for use
 * throughout this module. This circumvents endless calls to malloc.
 * Any function is free to use the buffer. We note that the buffer is
 * *not* assumed to have any default value. Those functions needing a
 * default value should remember to set it before using the buffer.  It
 * is perfectly legitimate for the user to resize the buffer when needed
 * but it is their responsibility to update the buffer size. Here
 * _BUFFSIZE records the size of the memory block pointed to by permbuff
 * in bytes (where byte means sizeof(char)). We note that the buffer
 * should never be shrunk, this simply just doesn't make sense. */
static size_t _BUFFSIZE;
static char *permbuff;

/* Utility function for resizing the buffer. */
static int
resize_permbuff(size_t newsize)
{
    char *tmp;

    if (newsize > _BUFFSIZE) {
        tmp = PyMem_Resize(permbuff, char, newsize);
        if (tmp == NULL) {
            PyErr_SetNone(PyExc_MemoryError);
            return -1;
        }
        permbuff = tmp;
        _BUFFSIZE = newsize;
    }

    return 0;
}

/* This macro, taken from the GAP3 source code, returns the image of I
 * under the permutation PT of degree DG. Don't use this with arguments
 * that have side effects! */
#define IMAGE(I, PT, DG) ( (I) < (DG) ? (PT)[(I)] : (I) )

/* This setup is fairly explanatory. We have array is a pointer to the
 * underlying data, degree describes the length of the list, and type is
 * either 1, 2 or 4 depending upon where the data in the array is stored
 * as either UINT1, UINT2 or UINT4. This is useful for switches. */
typedef struct {
    PyObject_HEAD
    char *array;
    UINT1 type;
    unsigned long degree;
} Perm;

static PyTypeObject Perm_Type; /* Forward */

/* These macros are useful for checking objects are Perms and for
 * accessing the members of perms when they have the PyObject type.  */
#define Perm_Check(v)  ((v)->ob_type == &Perm_Type)
#define Perm_array(v)  (((Perm *) (v))->array)
#define Perm_type(v)  (((Perm *) (v))->type)
#define Perm_degree(v)  (((Perm *) (v))->degree)

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


/***** INCLUDES OF INTEGER TEMPLATES MUST OCCUR BEFORE THIS POINT *****/

static PyObject *
Perm_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    PyObject *initial;
    size_t len;

    if (type == &Perm_Type && !_PyArg_NoKeywords("permutat.Perm()", kwds)) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O:Perm", &initial)) {
        return NULL;
    }

    /* Check we got a sequence or a permutation. */
    if (Perm_Check(initial)) {
        UINT1 typ = Perm_type(initial);

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

    if (len > MAX_4) {
        PyErr_SetString(PyExc_ValueError,
            "permutation exceeded maximum size");
        return NULL;
    }

    /* Now get the appropriate type of permutation. */
    if (len <= MAX_1) {
        return newperm_1(type, len, initial);
    }
    else if (len <= MAX_2) {
        return newperm_2(type, len, initial);
    }
    else {
        return newperm_4(type, len, initial);
    }
}

static void
Perm_dealloc(Perm* self)
{
    if (self->array != NULL) {
        PyMem_Del(self->array);
    }
    Py_TYPE(self)->tp_free((PyObject *) self);
}

static PyObject *
Perm_repr(Perm *self)
{
    UINT1 typ = self->type;

    switch (typ) {
        case 1:
            return Perm_repr_1(self);
        case 2:
            return Perm_repr_2(self);
        case 4:
            return Perm_repr_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_cycledecomp(Perm *self)
{
    UINT1 typ = self->type;

    switch (typ) {
        case 1:
            return Perm_cycledecomp_1(self);
        case 2:
            return Perm_cycledecomp_2(self);
        case 4:
            return Perm_cycledecomp_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_cycletype(Perm *self, PyObject *args)
{
    int flag = 0;
    UINT1 typ = self->type;

    if (!PyArg_ParseTuple(args, "|p:cycletype", &flag)) {
        return NULL;
    }

    if (flag) {
        switch (typ) {
            case 1:
                return Perm_cycletype_part_1(self);
            case 2:
                return Perm_cycletype_part_2(self);
            case 4:
                return Perm_cycletype_part_4(self);
        }
    }
    else {
        switch (typ) {
            case 1:
                return Perm_cycletype_1(self);
            case 2:
                return Perm_cycletype_2(self);
            case 4:
                return Perm_cycletype_4(self);
        }
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}


static PyObject *
Perm_order(Perm *self, void *closure)
{
    UINT1 typ = self->type;

    switch (typ) {
        case 1:
            return Perm_order_1(self);
        case 2:
            return Perm_order_2(self);
        case 4:
            return Perm_order_4(self);
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_cycle(PyObject *self, PyObject *args)
{
    PyObject *pt_py;
    unsigned long pt;

    if (!PyArg_ParseTuple(args, "O!:cycle", &PyLong_Type, &pt_py)) {
        return NULL;
    }

    pt = PyLong_AsUnsignedLong(pt_py);
    if (pt >= MAX_4 || (pt == (unsigned long) -1 && PyErr_Occurred())) {
        PyErr_Format(PyExc_ValueError,
            "permutations are only defined on integers in the range "
            "[0, ..., %lu]", MAX_4 - 1);
        return NULL;
    }

    if (Perm_degree(self) <= MAX_1) {
        return cycletuple_1((Perm *) self, pt);
    }
    else if (Perm_degree(self) <= MAX_2) {
        return cycletuple_2((Perm *) self, pt);
    }
    else {
        return cycletuple_4((Perm *) self, pt);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_cycles(PyObject *self, PyObject *args)
{
    if (Perm_degree(self) <= MAX_1) {
        return Perm_cycles_1((Perm *) self);
    }
    else if (Perm_degree(self) <= MAX_2) {
        return Perm_cycles_2((Perm *) self);
    }
    else {
        return Perm_cycles_4((Perm *) self);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_tolist(PyObject *self, PyObject *args)
{
    if (Perm_degree(self) <= MAX_1) {
        return Perm_tolist_1((Perm *) self);
    }
    else if (Perm_degree(self) <= MAX_2) {
        return Perm_tolist_2((Perm *) self);
    }
    else {
        return Perm_tolist_4((Perm *) self);
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

static PyObject *
Perm_onsets(Perm *self, PyObject *args)
{
    PyObject *list;
    Py_ssize_t len;
    size_t ulen;
    UINT1 typ = Perm_type(self);

    if (!PyArg_ParseTuple(args, "O!:onsets", &PyList_Type, &list)) {
        return NULL;
    }

    len = PyList_Size(list);
    if (len == -1) {
        return NULL;
    }

    ulen = (size_t) len;
    if (ulen > MAX_4) {
        PyErr_SetString(PyExc_ValueError,
            "The list contains too many sets to define a permutation.");
        return NULL;
    }

    if (ulen <= MAX_1) {
        switch (typ) {
            case 1:
                return Perm_onsets_11(self, list, ulen);
            case 2:
                return Perm_onsets_21(self, list, ulen);
            case 4:
                return Perm_onsets_41(self, list, ulen);
        }
    }
    else if (ulen <= MAX_2) {
        switch (typ) {
            case 1:
                return Perm_onsets_12(self, list, ulen);
            case 2:
                return Perm_onsets_22(self, list, ulen);
            case 4:
                return Perm_onsets_42(self, list, ulen);
        }
    }
    else {
        switch (typ) {
            case 1:
                return Perm_onsets_14(self, list, ulen);
            case 2:
                return Perm_onsets_24(self, list, ulen);
            case 4:
                return Perm_onsets_44(self, list, ulen);
        }
    }

    /* If we got here something went wrong. */
    PyErr_BadInternalCall();
    return NULL;
}


/* The following hash function is taken from the Python tuple hash
 * function. */
static Py_hash_t
permhash(Perm *self)
{
    unsigned long deg = self->degree;

    if (deg <= MAX_1) {
        return permhash_1(self);
    }
    else if (deg <= MAX_2) {
        return permhash_2(self);
    }
    else {
        return permhash_4(self);
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
Mod_cyclestoperm(PyObject *self, PyObject *args)
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
            0L, MAX_4 - 2);
        Py_DECREF(max_pylong);
        return NULL;
    }
    Py_DECREF(max_pylong);

    /* From this point there are no living Python objects. */
    if (max > MAX_4 - 2) {
        PyErr_Format(PyExc_ValueError,
            "elements of a cycle must lie in the range %lu to %lu",
            0L, MAX_4 - 2);
        return NULL;
    }

    if (max <= MAX_1) {
        return cyclestoperm_1(args, max + 1, numcycles);
    }
    else if (max <= MAX_2) {
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
Mod_restrictedperm(PyObject *self, PyObject *args)
{
    PyObject *perm, *list, *max_pylong = PyLong_FromLong(0L);
    Py_ssize_t i, len;
    unsigned long max, deg;
    UINT1 typ;

    if (max_pylong == NULL) {
        return NULL;
    }

    if (!PyArg_ParseTuple(args, "O!O:restrictedperm",
                         &Perm_Type, &perm, &list)) {
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
    if ((max == (unsigned long) -1 && PyErr_Occurred()) || max > MAX_4 - 2) {
        Py_DECREF(max_pylong);
        PyErr_Format(PyExc_ValueError,
            "The sequence must only contain integers between %lu and %lu.",
            0L, MAX_4 - 2);
        return NULL;
    }
    Py_DECREF(max_pylong);

    deg = max + 1;
    typ = Perm_type(perm);


    if (deg <= MAX_1) {
        switch (typ) {
            case 1:
                return restrictedperm_11(perm, list, len, deg);
            case 2:
                return restrictedperm_12(perm, list, len, deg);
            case 4:
                return restrictedperm_14(perm, list, len, deg);
        }
    }
    else if (deg <= MAX_2) {
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
                return PERM_CAT(PERM_CAT(Perm_, X), _11)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _12)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _14)(L, R);     \
        }                                                           \
    }                                                               \
    case 2: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _21)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _22)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _24)(L, R);     \
        }                                                           \
    }                                                               \
    case 4: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _41)(L, R);     \
            case 2:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _42)(L, R);     \
            case 4:                                                 \
                return PERM_CAT(PERM_CAT(Perm_, X), _44)(L, R);     \
        }                                                           \
    }                                                               \
}

static PyObject *
Perm_mul(PyObject *L, PyObject *R)
{
    UINT1 typL, typR;

    if (Perm_Check(L) && Perm_Check(R)) {
        typL = Perm_type(L);
        typR = Perm_type(R);

        _BIFUNC_RETURN(mul_perm_perm)
    }
    else if (Perm_Check(L) && PyTuple_Check(R)) {
        PyObject *max_pylong;
        unsigned long max;
        Py_ssize_t i, cyclen = PyTuple_Size(R);

        typL = Perm_type(L);
    
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
                0L, MAX_4 - 2);
            Py_DECREF(max_pylong);
            return NULL;
        }
        Py_DECREF(max_pylong);

        /* From this point there are no living Python objects. */
        if (max > MAX_4 - 2) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX_4 - 2);
            return NULL;
        }

        if (max <= MAX_1) {
            R = cycletoperm_1(R, max + 1, cyclen);
            typR = 1;
        }
        else if (max <= MAX_2) {
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
    else if (PyTuple_Check(L) && Perm_Check(R)) {
        PyObject *max_pylong;
        unsigned long max;
        Py_ssize_t i, cyclen = PyTuple_Size(L);

        typR = Perm_type(R);
    
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
                0L, MAX_4 - 2);
            Py_DECREF(max_pylong);
            return NULL;
        }
        Py_DECREF(max_pylong);

        /* From this point there are no living Python objects. */
        if (max > MAX_4 - 2) {
            PyErr_Format(PyExc_ValueError,
                "elements of a cycle must lie in the range %lu to %lu",
                0L, MAX_4 - 2);
            return NULL;
        }

        if (max <= MAX_1) {
            L = cycletoperm_1(L, max + 1, cyclen);
            typL = 1;
        }
        else if (max <= MAX_2) {
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
Perm_div(PyObject *L, PyObject *R)
{
    if (Perm_Check(L) && Perm_Check(R)) {
        UINT1 typL = Perm_type(L);
        UINT1 typR = Perm_type(R);

        _BIFUNC_RETURN(div_perm_perm)
    }
    else if (PyLong_Check(L) && Perm_Check(R)) {
        UINT1 typR = Perm_type(R);

        switch (typR) {
            case 1:
                return Perm_div_long_perm_1(L, R);
            case 2:
                return Perm_div_long_perm_2(L, R);
            case 4:
                return Perm_div_long_perm_4(L, R);
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
Perm_pow(PyObject *L, PyObject *R, PyObject *Z)
{
    UINT1 typL;

    if (!Perm_Check(L) || !PyLong_Check(R)) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    typL = Perm_type(L);

    switch (typL) {
        case 1:
            return Perm_pow_1(L, R);
        case 2:
            return Perm_pow_2(L, R);
        case 4:
            return Perm_pow_4(L, R);
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
Perm_xor(PyObject *L, PyObject *R)
{
    if (Perm_Check(L) && Perm_Check(R)) {
        UINT1 typL = Perm_type(L);
        UINT1 typR = Perm_type(R);

        _BIFUNC_RETURN(xor_perm_perm)
    }
    else if (PyLong_Check(L) && Perm_Check(R)) {
        UINT1 typR = Perm_type(R);

        switch (typR) {
            case 1:
                return Perm_xor_long_perm_1(L, R);
            case 2:
                return Perm_xor_long_perm_2(L, R);
            case 4:
                return Perm_xor_long_perm_4(L, R);
        }
    }
    else {
        Py_RETURN_NOTIMPLEMENTED;
    }

    /* If we got here something is wrong. */
    PyErr_BadInternalCall();
    return NULL;
}

/* Note that the variables c, L, R, typL, typR are not part of the macro
 * definition. It is assumed that they have already been defined. */
#define _BIFUNC_CMP(X)                                              \
switch (typL) {                                                     \
    case 1: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _11)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _12)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _14)(L, R);        \
        }                                                           \
        break;                                                      \
    }                                                               \
    case 2: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _21)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _22)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _24)(L, R);        \
        }                                                           \
        break;                                                      \
    }                                                               \
    case 4: {                                                       \
        switch (typR) {                                             \
            case 1:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _41)(L, R);        \
                break;                                              \
            case 2:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _42)(L, R);        \
                break;                                              \
            case 4:                                                 \
                c = PERM_CAT(PERM_CAT(Perm_, X), _44)(L, R);        \
        }                                                           \
    }                                                               \
}

static PyObject *
Perm_richcompare(PyObject *L, PyObject *R, int op)
{
    PyObject *res;
    int c;
    UINT1 typL, typR;

    if (!Perm_Check(L) || !Perm_Check(R)) {
        Py_RETURN_NOTIMPLEMENTED;
    }

    typL = Perm_type(L);
    typR = Perm_type(R);

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


static PyNumberMethods Perm_Type_AsNumber = {
     0,                         /* nb_add */
     0,                         /* nb_subtract */
     (binaryfunc) Perm_mul,     /* nb_multiply */
     0,                         /* nb_remainder */
     0,                         /* nb_divmod */
     (ternaryfunc) Perm_pow,    /* nb_power */
     0,                         /* nb_negative */
     0,                         /* nb_positive */
     0,                         /* nb_absolute */
     0,                         /* nb_bool */
     0,                         /* nb_invert */
     0,                         /* nb_lshift */
     0,                         /* nb_rshift */
     0,                         /* nb_and */
     (binaryfunc) Perm_xor,     /* nb_xor */
     0,                         /* nb_or */
     0,                         /* nb_int */
     0,                         /* *nb_reserved */
     0,                         /* nb_float */
     0,                         /* nb_inplace_add */
     0,                         /* nb_inplace_subtract */
     0,                         /* nb_inplace_multiply */
     0,                         /* nb_inplace_remainder */
     0,                         /* nb_inplace_power */
     0,                         /* nb_inplace_lshift */
     0,                         /* nb_inplace_rshift */
     0,                         /* nb_inplace_and */
     0,                         /* nb_inplace_xor */
     0,                         /* nb_inplace_or */
     (binaryfunc) Perm_div,     /* nb_floor_divide */
     (binaryfunc) Perm_div,     /* nb_true_divide */
     0,                         /* nb_inplace_floor_divide */
     0,                         /* nb_inplace_true_divide */
     0,                         /* nb_index */
     0,                         /* nb_matrix_multiply */
     0                          /* nb_inplace_matrix_multiply */
};


static PyMethodDef Perm_methods[] = {
    {"cycledecomp", (PyCFunction)Perm_cycledecomp, METH_NOARGS,
     "Give the decomposition of the permutation into disjoint cycles."},
    {"cycletype", (PyCFunction)Perm_cycletype, METH_VARARGS,
     "Give the cycle type of the permutation."},
    /* {"order", (PyCFunction)Perm_order, METH_NOARGS,
     "Compute the multiplicative order of the permutation."},*/
    {"cycle", (PyCFunction)Perm_cycle, METH_VARARGS,
     "Return the cycle of the permutation which moves the given point."},
    {"cycles", (PyCFunction)Perm_cycles, METH_NOARGS,
     "Return a list containing all the cycles of the permutation."},
    {"tolist", (PyCFunction)Perm_tolist, METH_NOARGS,
     "Return the permutation as a Python list."},
    {"onsets", (PyCFunction)Perm_onsets, METH_VARARGS,
     "Return the induced permutation acting on a list of sets."},
    {NULL} /* Sentinel */
};

static PyMemberDef Perm_members[] = {
    {"degree", T_ULONG, offsetof(Perm, degree), READONLY,
     "The minimum n such that the symmetric group S_n contains the permutation."},
    {NULL} /* Sentinel */
};

static PyGetSetDef Perm_getsetters[] = {
    {"order", (getter)Perm_order, NULL,
     "The multiplicative order of the permutation."},
    {NULL} /* Sentinel */
};


static PyMethodDef permutat_methods[] = {
    {"cyclestoperm", (PyCFunction)Mod_cyclestoperm, METH_VARARGS,
     "Return a permutation which is a product of the corresponding "
     "disjoint cycles."},
    {"restrictedperm", (PyCFunction)Mod_restrictedperm, METH_VARARGS,
     "Return the restriction of a permutation to a subset."},
    {NULL} /* Sentinel */
};


static PyTypeObject Perm_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "permutat.Perm",             /* tp_name */
    sizeof(Perm),                /* tp_basicsize */
    0,                           /* tp_itemsize */
    (destructor) Perm_dealloc,   /* tp_dealloc */
    0,                           /* tp_print */
    0,                           /* tp_getattr */
    0,                           /* tp_setattr */
    0,                           /* tp_reserved */
    (reprfunc) Perm_repr,        /* tp_repr */
    &Perm_Type_AsNumber,         /* tp_as_number */
    0,                           /* tp_as_sequence */
    0,                           /* tp_as_mapping */
    (hashfunc)permhash,          /* tp_hash  */
    0,                           /* tp_call */
    0,                           /* tp_str */
    0,                           /* tp_getattro */
    0,                           /* tp_setattro */
    0,                           /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT |
        Py_TPFLAGS_BASETYPE,     /* tp_flags */
    "Perm objects",              /* tp_doc */
    0,                           /* tp_traverse */
    0,                           /* tp_clear */
    Perm_richcompare,            /* tp_richcompare */
    0,                           /* tp_weaklistoffset */
    0,                           /* tp_iter */
    0,                           /* tp_iternext */
    Perm_methods,                /* tp_methods */
    Perm_members,                /* tp_members */
    Perm_getsetters,             /* tp_getset */
    0,                           /* tp_base */
    0,                           /* tp_dict */
    0,                           /* tp_descr_get */
    0,                           /* tp_descr_set */
    0,                           /* tp_dictoffset */
    0,                           /* tp_init */
    0,                           /* tp_alloc */
    Perm_new,                    /* tp_new */
};

static PyModuleDef permutat = {
    PyModuleDef_HEAD_INIT,
    "permutat",
    "Module giving python permutation arithmetic.",
    -1,
    permutat_methods, NULL, NULL, NULL, NULL
};

PyMODINIT_FUNC
PyInit_permutat(void)
{
    PyObject* m;

    if (PyType_Ready(&Perm_Type) < 0) {
        return NULL;
    }

    m = PyModule_Create(&permutat);
    if (m == NULL) {
        return NULL;
    }

    Py_INCREF(&Perm_Type);
    PyModule_AddObject(m, "Perm", (PyObject *)&Perm_Type);

    /* Initialise the buffer. */
    _BUFFSIZE = 1000 * sizeof(UINT2);
    permbuff = PyMem_New(char, _BUFFSIZE);
    if (permbuff == NULL) {
        return PyErr_NoMemory();
    }

    return m;
}
