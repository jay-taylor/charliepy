import numpy as np
import itertools
import sys
import collections

# define what is imported under import *
__all__ = ['cartanmat',
           'cartantotype']

def specialpairs(R):
    # For each positive root alpha we find the minimal simple root alpha_j
    # and positive root beta such that alpha = alpha_j + beta. Here R should be
    # a list consisting of the positive roots and n should be the number of
    # simple roots.
    #
    # Note that the ith entry of extspec is the extraspecial pair of the root
    # with position i+n.

    n = len(R.cartanmat)
    hts = [sum(r) for r in R.roots[:R.N]]
    mod = np.einsum("ij,jk,ik->i", R.roots[:R.N], R._Bform,
            R.roots[:R.N]).tolist()
    D = R._dictroots
    extspec = [None]*(R.N-n)
    pairs = dict()

    for i in range(n, R.N):
        first = True
        for j in range(i):
            tmp = tuple([a-b for a, b in zip(R.roots[i],R.roots[j])])
            if tmp in D:
                m = D[tmp]
                if j >= m:
                    if hts[j] > hts[m]:
                        break
                    else:
                        continue
                
                # Get the p and the q values for the root string. Here we use
                # the clasification of such root strings given in Lemma 4.1.1 of
                # [Car92].
                x, y = mod[m], mod[i]
                if x == y:
                    if (tuple([a-b for a,b in
                               zip(R.roots[m],R.roots[j])]) in D):
                        p, q = 1, 2
                    else:
                        p, q = 0, 1
                elif x == 2*y:
                    p, q = 0, 2
                elif 2*x == y:
                    p, q = 1, 1
                elif x == 3*y:
                    p, q = 0, 3
                elif 3*x == y:
                    p, q = 2, 1

                # Get the structure constant for the special pair.
                if first:
                    N = R.extraspecialsigns[i-n]*(p+1)
                    extspec[i-n] = [(j, m), p, q, N]
                    pairs.update({(j, m) : [i, p, q, N]})
                    first = False
                    continue
                else:
                    x, y = extspec[i-n][0][0], extspec[i-n][0][1]
                    tmp = tuple([a-b for a,b in
                        zip(R.roots[m],R.roots[x])])
                    if tmp in D:
                        k = D[tmp]
                        sgn = 1
                        if x > k:
                            p1 = (k,x)
                            sgn *= -1
                        else:
                            p1 = (x,k)
                        if j > k:
                            p2 = (k,j)
                            sgn *= -1
                        else:
                            p2 = (j,k)
                        t1 = sgn*mod[k]*mod[j]*pairs[p1][2]*pairs[p2][2]
                    else:
                        t1 = 0
                    
                    tmp = tuple([a-b for a,b in
                        zip(R.roots[j], R.roots[x])])
                    if tmp in D:
                        k = D[tmp]
                        sgn = 1
                        if x > k:
                            p1 = (k, x)
                            sgn *= -1
                        else:
                            p1 = (x, k)
                        if m > k:
                            p2 = (k, m)
                            sgn *= -1
                        else:
                            p2 = (m, k)
                        t2 = sgn*mod[k]*mod[m]*pairs[p1][2]*pairs[p2][2]
                    else:
                        t2 = 0
                    N = cmp(t1, t2)*(p+1)
                    N = p+1
                    pairs.update({(j, m) : [i, p, q, N]})
    
    R.specialpairs = pairs


def is_intcartanmat(C):
    """Determines whether an integral matrix is a Cartan matrix as defined in
    GP00.
    """
    n = len(C)
    vals = {0, 1, 2, 3, 4}
    for i in range(n):
        if C[i][i] != 2:
            return False
        for j in range(i+1, n):
            if C[i][j] > 0 or C[j][i] > 0:
                return False
            if C[i][j]*C[j][i] not in vals:
                return False
    return True

def cartanmat(*args, **kwargs):
    """
    If extras=True is passed as a keyword argument then the function
    returns a tuple containing: the Cartan matrix, the Cartan type of
    the matrix and the indices of the nonreduced parts of the matrix.

    Notes
    -----
    The Cartan matrices are as in [Bou02]. In particular, the Dynkin
    diagrams of the finite reduced crystallographic Cartan matrices are
    as follows:

         0   1   2        n-1          0   1            n-1
    An   o---o---o-- ... --o      Bn   o---o-- ... --o=>=o
    (n>0)                         (n>0)

                                                       o n-2
         0   1            n-1          0        n-4   /     
    Cn   o---o-- ... --o=<=o      Dn   o-- ... --o---o n-3  
    (n>0)                         (n>1)               \     
                                                       o n-1

         0   1           0   1   2   3          0   2   3   4   5
    G2   o-<-o       F4  o---o=>=o---o      E6  o---o---o---o---o
           6                                            |
                                                        o 1

         0   2   3   4   5   6           0   2   3   4   5   6   7
    E7   o---o---o---o---o---o       E8  o---o---o---o---o---o---o
                 |                               |
                 o 1                             o 1

    The Dynkin diagram of the unique non-reduced irreducible
    crystallographic Cartan matrix is given by:

                        0   1          n-1  n
                   BCn  o=>=o--- ... ---o=>=o

    The Dynkin diagrams of the affine Cartan matrices are:

          1   2       n-1  n           1 o
    A~n   o---o-- ...--o---o              \    3        n-1  n
    (n>1)  \      0       /      B~n     2 o---o-- ... --o=>=o
             -----o------        (n>2)    /
                                       0 o
    
                                       1 o                     o n-1
          0   1        n-1  n             \    3        n-3   /
    C~n   o=>=o-- ... --o=<=o    D~n     2 o---o-- ... --o---o n-2
    (n>1)                        (n>3)    /                   \ 
                                       0 o                     o n

          0   2   1           1   2   0           0   1   2   3   4
    B~2   o=>=o=<=o     G~2   o-<-o---o     F~4   o---o---o=>=o---o
                                6

          1   3   4   5   6          0   1   3   4   5   6   7
    E~6   o---o---o---o---o    E~7   o---o---o---o---o---o---o
                  |                              |
                  o 2                            o 2
                  |
                  o 0                1   3   4   5   6   7   8   0
                               E~8   o---o---o---o---o---o---o---o
                                             |
                                             o 2

    """
    rank = 0
    pairs = []
    ctyp = []
    try:
        extras = kwargs['extras']
    except KeyError:
        extras = False

    for typ in zip(*(iter(args),)*2):
        if typ[0][-1] == "~":
            nam = typ[0][:-1] + "t"
        else:
            nam = typ[0][:]
        try:
            rank += typ[1]
            if typ[0][-1] == "~" or typ[0] == "BC":
                rank += 1
                pairs.append((nam, typ[1]+1, False))
            else:
                pairs.append((nam, typ[1], False))
        except TypeError:
            rank += len(typ[1])
            pairs.append((nam, typ[1], True))

    cmat = np.zeros((rank, rank), dtype='int')

    if rank == 0:
        if extras:
            return cmat, [], []
        else:
            return cmat

    inds = collections.deque(range(rank))
    nonred = []
    for typ in pairs:
        sl = []
        if typ[2]: # Given the indices.
            for i in typ[1]:
                try:
                    inds.remove(i)
                except ValueError:
                    raise ValueError("Indices not aligned!")
                sl.append(i)
        else:
            for i in range(typ[1]):
                sl.append(inds.popleft())

        ctyp.append([typ[0],list(sl)])
        
        if sl:
            tmp = np.ix_(sl, sl)
            cmat[tmp] = sys.modules[
                    'charliepy.data.typ1' + typ[0]].cartanmat(len(sl))
            if typ[0] == 'BC':
                nonred.append([sl[0], sl[-1]])

    if inds:
        raise ValueError("Indices not aligned!")

    if extras:
        return cmat, ctyp, nonred 
    else:
        return cmat

def roots_finintcartan(C, nonred=[]):
    """
    Returns the roots of an integral Cartan matrix of finite type.

    Recall that every root is a Z-linear combination of the simple
    roots. Each root is thus represented by a tuple where the ith entry
    of the tuple gives the coefficient of the simple root alpha_i, with
    the indexing of the simple roots determined by the Cartan matrix.
    """
    # We use the following inductive algorithm to determine the roots.
    # This uses the assertions in Ch. VI, S1, no.  3, Proposition 9 of
    # [Bou02].
    #
    # Consider two roots alpha and beta such that alpha != beta.  The
    # set of integers i such that beta + i.alpha is a root is an
    # interval [-q, p] in Z containing 0. We call [-q, p] the
    # alpha-interval through beta. Clearly we have p-q is contained in
    # [-q, p] and moreover we have
    #
    #     p - q = -<beta, alpha*>.
    #
    # In particular, beta - <beta, alpha*>.alpha is a root. Now if beta
    # + alpha is a root then the alpha-interval through beta + alpha is
    # simply [-q-1, p-1]. As p - q - 1 is contained in this interval we
    # therefore have
    #
    #     beta - (<beta, alpha*> + 1).alpha
    #
    # is a root. Conversely if this is a root then we must have
    #
    #     -q <= p - q - 1  -->  1 <= p
    #
    # so beta + alpha is a root.
    #
    # This observation gives us an inductive approach to determining the
    # positive roots from the Cartan matrix. We proceed by induction on
    # the height of the root. The height 1 and 2 roots are easy to write
    # down. We then use the above characterisation to determine whether
    # adding a simple root gives a new positive root.
    #
    # When dealing with non-reduced root systems we always ensure that
    # each root has a unique expression by never using simple roots of
    # the form 2alpha_i.
    n = C.shape[0]
    droots = list(zip(*nonred))
    if not droots:
        droots = [[], []]
    inds = [i for i in range(n) if i not in droots[1]]
    I = np.identity(n, dtype='int8')

    # We start with the simple roots of height 1.
    roots = [[I[i].tolist() for i in inds]]

    # Now the roots of height 2.
    roots.append(
        [(I[i] + I[j]).tolist() for i in inds for j in range(i, n)
         if (j == i and i in droots[0])
         or (j != i and C[i][j] != 0 and j not in droots[1])]
    )

    # Now proceed by induction to obtain all roots.
    i = 1
    while roots[i]:
        roots.append([])
        prod = np.dot(roots[i], C).tolist()

        for j in range(len(roots[i])):
            r = roots[i][j]
            for k in inds:
                p = prod[j][k]

                # If p is negative then alpha+alpha_k is definitely a
                # root by Ch. VI, S1, no. 3, Corollary to Theorem 1 of
                # [Bou02].
                if p < 0:
                    r[k] += 1
                    if r not in roots[i+1]:
                        roots[i+1].append(list(r))
                    r[k] -= 1

                # If i <=p or r[k] <= p then we cannot have alpha -
                # (p+1)alpha_k is a root. After this all cheap tests are
                # done.
                elif i > p and r[k] > p:
                    r[k] -= p + 1
                    if r in roots[i-p-1]:
                        r[k] += p + 2
                        if r not in roots[i+1]:
                            roots[i+1].append(list(r))
                        r[k] -= 1
                    else:
                        r[k] += p + 1
        i += 1

    # We now just flatten the roots to a single list and add the
    # negative roots.  Note we have to create the tuples here as they
    # are immutable so we wouldn't be able to change their values in the
    # above loop.
    roots = [tuple(x) for x in itertools.chain.from_iterable(roots)]
    return tuple(roots + [tuple(-i for i in r) for r in roots])

def cartantotype(C):
    """The input is an integral Cartan matrix of finite type and the output is
    the type of the Cartan matrix. It is assumed that the input is genuinely a
    Cartan matrix.

    If the optional flags finite or crystal are set to True then the function
    will return False if the matrix is not finite or crystalographic
    respectively.
    """
    n = len(C)
    if not n:
        return [["A",[]]]

    # nonzero[i] gives all j such that C[i][j] != 0.
    nonzero = [np.nonzero(C[i])[0].tolist() for i in range(n)]
    blocks = []
    tmp = list(range(n))
    while tmp:
        orb = [tmp[0]]
        for i in orb:
            for j in tmp:
                if i in nonzero[j] and j in nonzero[i] and not j in orb:
                    orb.append(j)
        for i in orb: 
            tmp.remove(i)
        orb.sort()
        blocks.append(orb)

    typ = []
    for b in blocks:
        tnonzero = [[b.index(j) for j in nonzero[i]] for i in b]
        ntyp = cartantotype_indec(C[np.ix_(b,b)], tnonzero)
        ntyp[1] = [b[i] for i in ntyp[1]]
        typ.append(ntyp)

    return typ


def cartantotype_indec(C, nonzero):
    """Determines the Cartan type when the matrix is indecomposable. In
    particular the corresponding Dynkin diagram of C is connected.
    """
    n = len(C)
    if n == 1: # A1
        return ["A", [0]]

    if n == 2: # Rank 2.
        if C[0][1] == -4 and C[1][0] == -1: return ["BC", [0,1]]
        if C[1][0] == -4 and C[0][1] == -1: return ["BC", [1,0]]
        if C[1][0] == -3 and C[0][1] == -1: return ["G", [0,1]]
        if C[0][1] == -3 and C[1][0] == -1: return ["G", [1,0]]
        if C[0][1] == -2 and C[1][0] == -1: return ["B", [0,1]]
        if C[1][0] == -2 and C[0][1] == -1: return ["B", [1,0]]
        if C[0][1] == -1 and C[1][0] == -1: return ["A", [0,1]]
        if C[0][1] == -2 and C[1][0] == -2: return ["A~", [0,1]]
        return ["U", [0,1]]
    
    # Must now have rank at least 3.
    ends = []
    bpoints = []
    sbpoints = []

    uniq = np.unique(C)
    sbond = np.array_equal(uniq, [-1, 0, 2]) # single bond
    dbond = np.array_equal(uniq, [-2, -1, 0, 2]) # double bond
    tbond = np.array_equal(uniq, [-3, -1, 0, 2]) # triple bond

    for i in range(n):
        if len(nonzero[i]) == 2: # End point
            ends.append(i)
        elif len(nonzero[i]) == 4: # Branch point
            bpoints.append(i)
        elif len(nonzero[i]) == 5: # D4~ branch point
            sbpoints.append(i)
        elif len(nonzero[i]) != 3: # Don't know what it is.
            return ["U", list(range(n))]

    # All vertices now have degree 2, 3, 4 or 5.
    if len(sbpoints) == 1 and sbond: # D4~
        if len(ends) == 4 and n == 5:
            typinds = sorted(ends)
            typinds.insert(3, sbpoints[0])
            return ["D~", typinds]
    if sbpoints:
        return ["U", list(range(n))]

    # All vertices now have degree 2, 3, or 4.
    if len(ends) == 4 and len(bpoints) == 2 and sbond:
        # Could be Dn~
        grpd = []
        for p in bpoints:
            grpd.append([j for j in nonzero[p] if j in ends])
        if (len(grpd[0]),len(grpd[1])) != (2,2):
            return ["U", list(range(n))]
        # Must now be Dn~
        a, b = bpoints[0], bpoints[1]
        grpd[0].sort()
        grpd[1].sort()
        typinds = grpd[0][:]
        typinds.append(a)
        i = [j for j in nonzero[a] if j != a and j not in grpd[0]][0]
        while i != b:
            typinds.append(i)
            i = [j for j in nonzero[i] if j not in typinds][0]
        typinds.append(b)
        return ["D~", [grpd[1][0]] + typinds + [grpd[1][1]]]
    
    if len(ends) == 3 and len(bpoints) == 1:
        # Note, we must now have precisely 3 straight line branches.
        a = bpoints[0]
        branches = []
        tinds = set([a])
        for e in ends:
            tinds.add(e)
            branch = [e]
            i = e
            while a not in nonzero[i]:
                i = [j for j in nonzero[i] if j not in tinds][0]
                tinds.add(i)
                branch.append(i)
            branches.append(branch)
        branches.sort(key=(lambda x: len(x)))
        blens = [len(branch) for branch in branches]

        if sbond: # En, En~, Dn.
            if blens == [1, 2, 2]: # E6
                if branches[1][0] > branches[2][0]:
                    branches[2], branches[1] = branches[1], branches[2]
                branches[2].reverse()
                branches[1].insert(1, branches[0][0])
                typinds = branches[1] + [a] + branches[2]
                return ["E", typinds]
            if blens == [1, 2, 3] or blens == [1, 2, 4]: # E7, E8
                branches[2].reverse()
                branches[1].insert(1, branches[0][0])
                return ["E", branches[1] + [a] + branches[2]]
            if blens == [2, 2, 2]: # E6~
                branches.sort(key=(lambda x: x[0]))
                typinds = [None]*n
                typinds[:3:2] = branches[0]
                typinds[1:4:2] = branches[1]
                typinds[4] = a
                typinds[5:] = branches[2][::-1]
                return ["E~", typinds]
            if blens == [1, 3, 3]: # E7~
                if branches[1][0] > branches[2][0]:
                    branches[2], branches[1] = branches[1], branches[2]
                typinds = [None]*n
                typinds[:2] = branches[1][:2]
                typinds[2] = branches[0][0]
                typinds[3] = branches[1][2]
                typinds[4] = a
                typinds[5:] = branches[2][::-1]
                return ["E~", typinds]
            if blens == [1, 2, 5]: # E8~
                typinds = [None]*n
                typinds[0] = branches[2][0]
                typinds[1] = branches[1][0]
                typinds[2] = branches[0][0]
                typinds[3] = branches[1][1]
                typinds[4] = a
                typinds[5:] = branches[2][-1:0:-1]
                return ["E~", typinds]
            if blens == [1, 1, 1]: # D4
                typinds = branches[0] + branches[1] + branches[2]
                typinds.sort()
                typinds.insert(2, a)
                return ["D", typinds]
            if blens[:2] == [1, 1]: # Dn
                typinds = branches[0] + branches[1]
                typinds.sort()
                return ["D", typinds + [a] + branches[2][::-1]]
            return ["U", list(range(n))]

        if dbond and blens == [1, 1, 1]: # B3~
            dbend = []
            for e in ends:
                if -2 in C[e]:
                    dbend.append(e)
            if len(dbend) != 1 or -2 in C[a]:
                return ["U", list(range(n))]
            e = dbend[0]
            ends.remove(e)
            ends.sort()
            typinds = [ends[0]] + [e, a] + [ends[1]]
            return ["B~", typinds]

        if dbond and blens[:2] == [1, 1]: # Bn~
            if np.where(C == -2) == (branches[2][0], branches[2][1]):
                if branches[0][0] > branches[1][0]:
                    branches[1], branches[0] = branches[0], branches[1]
                return ["B~", branches[0] + branches[2] + [a] + branches[1]]
            return ["U", list(range(n))]

    if bpoints:
        return ["U", list(range(n))]

    # All vertices now have degree 2 or 3.
    if len(ends) == 2: # An, Bn, Cn, BCn, F4, B2~, Cn~, F4~, G2~.
        # Note, we must now have a straight line.
        ends.sort()
        e = ends[0] # line starts as lex smallest end node.
        line = [e]
        while e != ends[1]:
            e = [j for j in nonzero[e] if j not in line][0]
            line.append(e)

        if sbond: # An.
            return ["A", line]
        if dbond:
            inds = zip(*np.where(C == -2))
            if len(inds) == 1: # Bn, Cn, F4, F4~.
                bond = inds[0]
                if n == 4 and bond == (line[2],line[1]):
                    return ["F", line]
                if n == 4 and bond == (line[1],line[2]):
                    return ["F", line[::-1]]
                if n == 5 and bond == (line[3],line[2]):
                    return ["F~", line]
                if n == 5 and bond == (line[1],line[2]):
                    return ["F~", line[::-1]]
                if bond == (line[0],line[1]):
                    return ["B", line]
                if bond == (line[n-1], line[n-2]):
                    return ["B", line[::-1]]
                if bond == (line[1],line[0]):
                    return ["C", line]
                if bond == (line[n-2],line[n-1]):
                    return ["C", line[::-1]]
            if len(inds) == 2: # B2~, Cn~.
                S = set(inds)
                if n == 3 and S == set([(line[1],line[0]), (line[1],line[2])]):
                        return ["B~", line]
                if n > 3 and S == set([(line[1], line[0]),
                                       (line[n-2], line[n-1])]):
                    return ["C~", [line[0]] + line[-1:0:-1]]
                if S == set([(line[0],line[1]), (line[1],line[2])]):
                    return ["BC", line]
                if S == set([(line[2],line[1]), (line[1],line[0])]):
                    return ["BC", line[::-1]]
        if tbond: # G2~.
            inds = zip(*np.where( C == -3))
            bond = inds[0]
            if n == 3 and len(inds) == 1:
                if bond == (line[2], line[1]):
                    return ["G~", line]
                if bond == (line[0], line[1]):
                    return ["G~", line[::-1]]

    if ends:
        return ["U", list(range(n))]

    # All vertices now have degree 3. Must be a circuit.
    i = 0
    typinds = [0]
    while len(typinds) < n:
        i = [j for j in nonzero[i] if j not in typinds][0]
        typinds.append(i)
    return ["A~", typinds]

