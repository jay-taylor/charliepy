#ifndef CLP_PERMMODULE_H
#define CLP_PERMMODULE_H
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Our permutations will be defined on at most 2**32-2 points. This
 * means that the degree of the permutation is at most 2**32-1, which is
 * the maximum value for an unsigned 32-bit integer. To save memory we
 * consider three different types of permutations, depending upon
 * whether the permutation acts on 2**8-1, 2**16-1, or 2**32-1 points.
 * For this we introduce the following types. */
typedef uint_least8_t ClpPerm_UINT1;
typedef uint_least16_t ClpPerm_UINT2;
typedef uint_least32_t ClpPerm_UINT4;

/* This setup is fairly explanatory. We have array is a pointer to the
 * underlying data, degree describes the length of the list, and rank is
 * either 1, 2 or 4 depending upon whether the data in the array is stored
 * as either UINT1, UINT2 or UINT4. This is useful for switches. */
typedef struct {
    PyObject_HEAD
    char *array;
    unsigned short rank;
    unsigned long degree;
} ClpPermObject;

/* These macros are useful for checking objects are ClpPerms and for
 * accessing the members of perms when they have the PyObject type.  */
#define ClpPerm_Check(ob)  (Py_TYPE(ob) == &ClpPerm_Type)
#define ClpPerm_ARRAY(ob)  (((ClpPermObject *) (ob))->array)
#define ClpPerm_RANK(ob)  (((ClpPermObject *) (ob))->rank)
#define ClpPerm_DEGREE(ob)  (((ClpPermObject *) (ob))->degree)

/* Public C API table */
typedef struct {
    char *(*buffer)(void);
    int (*resizebuffer)(size_t);
    void (*zerobuffer)(unsigned long);
    PyObject *(*newperm)(unsigned long);
    PyObject *(*identity)(unsigned long);
    PyObject *(*cycledecomp)(ClpPermObject *);
    PyObject *(*cycle)(PyObject *, PyObject *);
    PyObject *(*order)(ClpPermObject *, void *);
    PyObject *(*cycletype)(ClpPermObject *, void *);
    PyObject *(*cyclestructure)(ClpPermObject *, void *);
    PyObject *(*cycles)(PyObject *, PyObject *);
    PyObject *(*tolist)(PyObject *, PyObject *);
    PyObject *(*onsets)(ClpPermObject *, PyObject *);
    PyObject *(*permute)(ClpPermObject *, PyObject *);
    PyObject *(*restrictedperm)(PyObject *, PyObject *);
} _ClpPermAPIMethods;

#ifndef CLP_PERM_MODULE
/* Use this section in client modules. */
static _ClpPermAPIMethods *_ClpPerm_API;

/* Import the API table from permutat.
 *
 * Note: We can't use PyCapsule_Import here because the permutat module
 * is a submodule of the charliepy package. This is too complicated for
 * PyCapsule_Import so we first have to import the module and then get
 * the C API attribute directly. */
static int
import_permutat(void)
{
    PyObject *perm_mod, *api_obj;

    /* Import the permutation module. Note this returns permutat, not
     * charliepy, because of the '.' in the name. */
    perm_mod = PyImport_ImportModule("charliepy.permutat");
    if (perm_mod == NULL) {
        return -1;
    }

    /* Get the API capsule object. */
    api_obj = PyObject_GetAttrString(perm_mod, "_C_API");
    Py_DECREF(perm_mod);
    if (api_obj == NULL) {
        return -1;
    }

    _ClpPerm_API = 
        (_ClpPermAPIMethods *) PyCapsule_GetPointer(api_obj, "permutat._C_API");
    return (_ClpPerm_API != NULL) ? 0 : -1;
}

/* Macros to implement the API. */
#define ClpPerm_Buffer() (_ClpPerm_API->buffer)()
#define ClpPerm_ResizeBuffer(x) (_ClpPerm_API->resizebuffer)(x)
#define ClpPerm_ZeroBuffer(x) (_ClpPerm_API->zerobuffer)(x)
#define ClpPerm_New(x) (_ClpPerm_API->newperm)(x)
#define ClpPerm_Identity(x) (_ClpPerm_API->identity)(x)
#define ClpPerm_CycleDecomp(x) (_ClpPerm_API->cycledecomp)(x)
#define ClpPerm_Cycle(x, y) (_ClpPerm_API->cycledecomp)(x, y)
#define ClpPerm_Order(x, y) (_ClpPerm_API->order)(x, y)
#define ClpPerm_CycleType(x, y) (_ClpPerm_API->cycletype)(x, y)
#define ClpPerm_CycleStructure(x, y) (_ClpPerm_API->cyclestructure)(x, y)
#define ClpPerm_Cycles(x, y) (_ClpPerm_API->cycles)(x, y)
#define ClpPerm_ToList(x, y) (_ClpPerm_API->tolist)(x, y)
#define ClpPerm_OnSets(x, y) (_ClpPerm_API->onsets)(x, y)
#define ClpPerm_Permute(x, y) (_ClpPerm_API->permute)(x, y)
#define ClpPerm_RestrictedPerm(x, y) (_ClpPerm_API->restrictedperm)(x, y)
#endif

#ifdef __cplusplus
}
#endif
#endif
