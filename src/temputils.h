#ifndef _TEMPUTILS_H
#define _TEMPUTILS_H

#define _PERM_CAT(A, B) A##B
#define PERM_CAT(A, B) _PERM_CAT(A, B)
#define PERM_UCAT(A, B) PERM_CAT(A, PERM_CAT(_, B))

#endif
