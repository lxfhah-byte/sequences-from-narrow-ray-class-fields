Sage code: (even dimension cases)
\begin{lstlisting}
q=11
F=GF(q)
E=EllipticCurve(F,[0,0,0,1,1])
Pts=[p for p in E if not p.is_zero()]
D=Pts
n=len(D)
single_k = 4

def eval_basis(k,ell,eta,pt):
    xp,yp = pt[0], pt[1]
    mi=k-2//2
    mj=(k-4)//2
    vec=[]
    for i in range(mi+1):
        vec.append(xp^i)
    for j in range(mj+1):
        vec.append((xp^j)*yp)
    e=k//2
    vec.append(xp^e+eta*(xp^mi)*yp)
    return vec

def codeGen(k,ell,eta):
    rows=[]
    for p in D:
        rows.append(eval_basis(k,ell,eta,p))
    M=matrix(F,rows)
    Mr=M.rref()
    new_rows=[]
    for r in Mr:
        if not r.is_zero():
            new_rows.append(r)
    if len(new_rows)==0:
        return None
    Mclean=matrix(F,new_rows)
    return LinearCode(Mclean)

cnt_k = 0
if single_k +1 < n:
    Lmax=single_k//2
    for l in range(Lmax+1):
        for t in F:
            if t == 0:
                continue
            C=codeGen(single_k,l,t)
            if not C:
                continue
            d=C.minimum_distance()
            nd,kc=C.length(),C.dimension()
            sd=nd-kc+1
            if d == sd:
                cnt_k +=1
print(cnt_k,n)
\end{lstlisting}
(odd dimension cases)
\begin{lstlisting}
q=5
    F = GF(q)
E = EllipticCurve(F, [0,0,0,1,1])
Pts = [p for p in E if not p.is_zero()]
D = Pts
n = len(D)
single_k = 3

def eval_basis(k,ell,eta,pt):
    xp,yp = pt[0], pt[1]
    mi = (k-1)//2
    mj = (k-5)//2
    vec = []
    for i in range(mi+1):
        vec.append(xp^i)
    for j in range(mj+1):
        vec.append((xp^j)*yp)
    e = (k+1)//2
    s=(k-3)//2
    vec.append( (xp^s)*yp+eta*(xp^e))
    return vec

def codeGen(k,ell,eta):
    rows = []
    for p in D:
        rows.append(eval_basis(k,ell,eta,p))
    M = matrix(F, rows)
    Mr = M.rref()
    new_rows = []
    for r in Mr:
        if not r.is_zero():
            new_rows.append(r)
    if len(new_rows) == 0:
        return None
    Mclean = matrix(F, new_rows)
    return LinearCode(Mclean)

cnt_k = 0
if single_k + 1 < n:
    Lmax = single_k//2
    for l in range(Lmax+1):
        for t in F:
            if t == 0:
                continue
            C = codeGen(single_k,l,t)
            if not C:
                continue
            d = C.minimum_distance()
            nd,kc = C.length(), C.dimension()
            sd = nd - kc + 1
            if d == sd:
                cnt_k += 1
print(cnt_k,n)
