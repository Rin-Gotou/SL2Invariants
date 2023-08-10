def prod(pollist,w):
    f = 1
    for i in range(len(w)):
        if w[i] > 0:
            f = f * pollist[i]^w[i]
    return f

def picker(pollist,kermat):
    l = len(pollist)
    m = len(kermat)
    ret = []
    idxs = []
    r = 0
    for c in range(l):
        if r >= m or kermat[r][c] == 0:
            ret.append(pollist[c])
            idxs.append(c)
        else:
            r = r+1
    print(idxs)
    return idxs,ret

'''
decompR22(f):
input:
    f: invariant in k[V_4+V_2+V_0+V_1^o]

output:
    multi-degree decomposition of f in dict form
    keys: the tuple of multi-degree in f
    value: each multi-degree part of f
'''

def decompR22(f):
    ret = dict()
    sepmat = matrix([[1,0,0],[1,0,0],[1,0,0],[1,0,0],[1,0,0],[0,1,0],[0,1,0],[0,1,0],[0,0,1],[0,0,0]])
    for t in f:
        c,m = t
        dbig = m.degrees()
        dpart = tuple(dbig * sepmat)
        if not dpart in ret.keys():
            ret[dpart] = 0
        ret[dpart] = ret[dpart] + c*m
    return ret
'''
lindec(invlist, f, degvecs, deg, expvals, omitmons)
linrels(invlist, degvecs, deg, expvals, omitmons)

input:
    invlist = [a_1,a_2,...,a_n]: the list of generating invariants
    f: an invariant
    degvec = (w_1,...,w_n): multidegrees of invariants
    deg: degree of f
    expvals = [x_1,x_2,...,x_n]: variables to expression
    omitmons = [m_j]: list of monomials of x_i's, which are leading monomials of f_j(x_i) s.t. f_j(a_i) = 0.

output:
    lindec: g(x_1,...,x_n) s.t. g(a_1,...,a_n) = f and any term aren't divisible by any m_j.
    linrels: g_j(x_1,...,x_n) s.t. g_j(a_1,...,a_n) = 0

method: by taking kernel of coefficient matrix (may not be completely efficient)
'''

def lindec(invlist, f, degvecs, deg, expvars, omitmons):
    mvecs = multdegvecs(deg,degvecs)
    expmons = [prod(expvars,v) for v in mvecs]
    if omitmons == []:
        imvecs = mvecs
    else:
        I = ideal(omitmons)
        imvecs = []
        for i in range(len(expmons)):
            if I.reduce(expmons[i]) != 0:
                imvecs.append(mvecs[i])
    expmons = [prod(expvars,v) for v in imvecs]
    pl = [f] + [prod(invlist,v) for v in imvecs]
    coes,mons = Sequence(pl).coefficient_matrix()
    kmat = kernel(coes).basis()
    vec = kmat[0][1:]
    ret = -vec * vector(expmons)
    return ret

def linrels(invlist, degvecs, deg, expvars, omitmons):
    mvecs = multdegvecs(deg,degvecs)
    expmons = [prod(expvars,v) for v in mvecs]
    if omitmons == []:
        imvecs = mvecs
    else:
        I = ideal(omitmons)
        imvecs = []
        for i in range(len(expmons)):
            if I.reduce(expmons[i]) != 0:
                imvecs.append(mvecs[i])
    expmons = [prod(expvars,v) for v in imvecs]
    pl = [prod(invlist,v) for v in imvecs]
    kmat = []
    if imvecs != []:
        coes,mons = Sequence(pl).coefficient_matrix()
        kmat = kernel(coes).basis()
    if kmat == []:
        return []
    ret = matrix(kmat) * vector(expmons)
    return list(ret)

