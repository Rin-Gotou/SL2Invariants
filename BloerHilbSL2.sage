
"""
hilbmat(S):
hilbfcn(S):
    input: S = [s_1,s_2,...]: sequence of nonnegative integers.
    Make the Bloer's matrix of the SL_2 representation V_S = V_s_1 + V_s_2 + ... .
    output:
        hilbmat(S): Bloer's matrix
        hilbfcn(S): multi-graded Hilbert series of V_S.
    Remark:
        Bloer's matrix for the SL_2 (root system A_1 and the positive root a_1):
            F_- := prod (1 - x^j t_i) for j: wt of V_s_i and (j,a) < 0
            F_b := B(b * F) for b >= 0
                B : Demazure operator, in this case B(x^i) := x^(abs(i+1)-1) , B(x^-1) := 0, QQ[t_i]-linear
            M := (coeff of x^a in F_b)_{a,b >= 0} is the Bloer's matrix.
        In program, 't0' is used instead of 'x' in this comment.
    Example:
        S = [3,2]
        V_3 : wt : {-3,-1,1,3}, V_2 : wt : {-2,0,2}
        F_- := (1-x^-3*t1) * (1-x^-1*t1) * (1-x^-2*t2)
"""

def hilbmat(S):
    n = len(S)
    ts = ["t"+str(i) for i in range(n+1)]
    RS = PolynomialRing(QQ,ts)
    deg = sum([ floor((s+1)/2) * floor((s+2)/2) for s in S])
    f = product([(1- RS.gen(i+1)*RS.gen(0)^j) for i in range(n) for j in range(S[i],0,-2)])
    fcoe = [f.coefficient({RS.gen(0):i}) for i in range(deg+1)]
    fcoe = [0 for i in range(deg)] + fcoe + [0 for i in range(deg+2)]
    Amat = matrix([[ fcoe[j-i+deg] - fcoe[i+j+deg+2] for i in range(deg+1)] for j in range(deg+1)])
    return Amat

def hilbfcn(S):
    n = len(S)
    ts = ["t"+str(i) for i in range(n+1)]
    RS = PolynomialRing(QQ,ts,order = "negdeglex(" + str(n) + "),neglex(1)")
    deg = sum([ floor((s+1)/2) * floor((s+2)/2) for s in S])
    f = product([(1- RS.gen(i+1)*RS.gen(0)^j) for i in range(n) for j in range(S[i],0,-2)])
    fcoe = [f.coefficient({RS.gen(0):i}) for i in range(deg+1)]
    Amat = matrix([[ (fcoe[j-i] if j-i >= 0 else 0) - (fcoe[i+j+2] if i+j+2 <= deg else 0) for i in range(deg+1)] for j in range(deg+1)])
    return factor(Amat.submatrix(1,1,deg,deg).det())/factor(Amat.det()* product([1 - RS.gen(i+1) if S[i] % 2 == 0 else 1 for i in range(n)]))


"""
coeffs(f,x): unused shorthand to make coefficient-list in one chosen variable x.
    Input:
        f : polynomial
        x : generator of f.par()
    Output:
        list [f0,f1,...] s.t. f = f0 + f1*x + f2*x^2 + ... .
There is an alternative way f.polynomial(x).coefficients(), but this makes an instance of a new polynomial ring.
"""
def coeffs(f,x):
    fdeg = f.degree(x)
    flist = [f.coefficient({x:i}) for i in range(fdeg+1)]
    return flist
