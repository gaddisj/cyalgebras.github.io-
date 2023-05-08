var('t')

#generic functions for matrices

def spec_rad(M):
	E=M.eigenvalues()
	return(max([abs(e) for e in E]))

def is_strongly_connected(M):
	d=len(M.rows())
	W=identity_matrix(d)
	z=zero_vector(d)
	e=(W.eigenvectors_right())[0][1]
	for i in range(d):
		v=M.row(i)-M[i,i]*e[i]
		if v==z:
			return(False)
		v=M.column(i)-M[i,i]*e[i]
		if v==z:
			return(False)
	return(True)

def is_normal(M):
	N=M.transpose()
	Y=M*N-N*M
	return(Y==0)

#CY specific functions

def mpoly(M,P,s):
	return((identity_matrix(M.nrows())-M*t+P*(M.transpose())*t^(s-1)-P*t^s).determinant())

def hilb_series(M,P,n,s):
	N=M.transpose()
	if s==3:
		H=[M^0,M,M^2-P*N]
		for i in range(3,n):
			H.append(M*H[i-1]-P*N*H[i-2]+P*H[i-3])
	if s==4:
		H=[M^0,M,M^2,N^3-P*N]
		for i in range(4,n):
			H.append(M*H[i-1]-P*N*H[i-3]+P*H[i-4])
	return(H)

def trunc_hilb_series(M,P,n):
	N=M.transpose()
	H=[M^0,M,M^2-P*N]
	for i in range(3,n):
		H.append(M*H[i-1]-P*N*H[i-2])
	return(H)

def new_cyc_list(deg):
	if deg==2:
		return([t^2-2*t+1,t^2+2*t+1,t^2+t+1,t^2-t+1,t^2+1])
	#since phi(n)>= sqrt(n/2), then we have the following upper bound for our polynomials
	dmax=2*(deg+1)^2
	list_of_polys=[]
	for i in range(2,dmax):
		cyp=cyclotomic_polynomial(i,t) 
		if 1 < cyp.degree() == deg:
			list_of_polys.append(cyp)
	return(list_of_polys)

def cyc_prods(L):
	l = len(L)
	if l == 1:
		return(new_cyc_list(2*L[0]))
	fac1 = new_cyc_list(2*L[0])
	if l == 2:
		fac2 = new_cyc_list(2*L[1])
		#prods = [a*b for a in fac1 for b in fac2]
	else:
		fac2 = cyc_prods(L[1:l])
	prods = [expand(a*b) for a in fac1 for b in fac2] 
	return(sorted(set(prods)))


def get_cycs(deg):
	list_of_prods=new_cyc_list(deg)
	parts=Partitions(Integer(deg/2)).list()
	del parts[0]
	for p in parts:
		list_of_prods=list_of_prods+cyc_prods(p)
	return(sorted(set(list_of_prods)))

def sort_cycs(deg):
	cyc_list=get_cycs(deg)
	srt=[(i,[]) for i in range(-deg,deg+1)]
	for p in cyc_list:
		cf=p.list()[1]
		srt[cf+deg][1].append(p)	
	return(srt)

def lam_polys(lam,P,s):
	det=P.determinant()
	if s==3: 
		if det==1:
			S=sort_cycs(8)[lam+12][1]
			new_list=[expand((1-t)^4*(poly)) for poly in S]
		if det==-1:
			lam=-lam
			S=sort_cycs(8)[lam+10][1]
			new_list=[expand((1+t)*(1-t)^3*(poly)) for poly in S]
	if s==4: 
		if det==1:
			S=sort_cycs(12)[lam+16][1]
			new_list=[expand((1-t)^4*(poly)) for poly in S]
		if det==-1:
			lam=-lam
			S=sort_cycs(12)[lam+14][1]
			new_list=[expand((1+t)*(1-t)^3*(poly)) for poly in S]
	return(new_list)

def uniq(lst):
    last = object()
    for item in lst:
        if item == last:
            continue
        yield item
        last = item

def sort_and_deduplicate(l):
    return list(uniq(sorted(l, reverse=True)))


def clear_dupes(mlist):
	G=SymmetricGroup(4)
	gmats=[g.matrix() for g in G]
	new_list=[]
	for N in mlist:
		Nperms=[Q*N*Q.inverse() for Q in gmats]
		intersect=[K for K in Nperms if K in new_list]
		if len(intersect)==0:
			new_list.append(N)
	return(new_list)

#Returns the mutation (in the sense of Fomin-Zelevinsky
#of the quiver with adjacency matrix M at vertex v
#note vertices are numbered 0,...,dim-1
def mutate(M,v):
	dim = len(M.rows())
	if dim!=len(M.columns()):
		return(print("error: not square"))
	if (v < 0) or (dim-1 < v): 
		return(print("error: invalid vertex"))
	for i in range(dim):
		if M[i,i]!=0:
			return(print("error: not loopless")) 
		for j in range(1,dim):
			if M[i, j]!=0 and M[j, i]!=0:
				return(print("error: two-cycles"))
	#Find arrows going into or out of vertex v, adds an arrows to M creating matrix K
	N=Matrix(dim)
	for i in range(dim):
		for j in range(dim):
			N[i, j]=M[i, v]*M[v, j]
	K=M + N
	#switches arrows in or out of vertex v
	U=M.column(v)
	V=M.row(v) 
	K[v,:]=U
	K[:,v]=V
	#remove 2-cycles
	for i in range(dim):
		for j in range(i+1,dim):
			a=K[i, j]
			b=K[j, i]
			if a!=0 and b!=0:
				if b < a: 
					K[i, j]=a - b
					K[j, i]= 0
				else: 
					K[i, j]=0
					K[j, i]=b - a
	return(K)