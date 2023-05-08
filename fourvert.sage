load('cypack.sage')

#These routines are specific to the 4-vertex quiver analysis.
#See paper:
#"Four-vertex quivers supporting twisted graded Calabi-Yau algebras" by Gaddis, Lamkin, Nguyen, Wright

###############################################################
#Four-cycle
###############################################################

def four_chk(s):
	P=Matrix([[0,1,0,0],[0,0,1,0],[0,0,0,1],[1,0,0,0]]) 
	mlist=[]
	d=6-s
	L=[(i,j,k,l) for i in range(d+1) for j in range(d+1) for k in range(d+1) for l in range(d+1) if i+j+k+l==d]
	for f in L:
		M=matrix.circulant(f) 
		q=mpoly(M,P,s)
		c=q.coefficient(t^(4*s-1))
		if q in lam_polys(c,P,s):
			mlist.append(matrix.circulant(f))
	return(clear_dupes(mlist))

###############################################################
#Three-cycle
###############################################################

def three_check():
	P=Matrix([[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]])
	cy1= -r - 3*y - 3
	cy2= -3*u*v + 3*r*y + 4*r + 12*y + 3
	L=[(a,b) for a in range(4) for b in range(4)]
	mat_list=[]
	for tuple in L:
		lam=cy1.subs(r==tuple[0],y==tuple[1])
		polys=lam_polys(lam,P,3)
		for poly in polys:
			beta=cy2.subs(r==tuple[0],y==tuple[1])
			cy3=(beta-poly.coefficient(t^10))/(-3)
			print(tuple,cy3)

def uv_evals(s):
	var('w x y u r')
	M=Matrix([[w,x,y,u],[y,w,x,u],[x,y,w,u],[u,u,u,r]]) 
	if s==3:
		L=[(a,b,c,d,e) for a in range(4) for b in range(5) for c in range(7) for d in range(7) for e in range(2) if a+b+c+d <= 6]
	if s==4:
		L=[(a,b,c,d,e) for a in range(3) for b in range(5) for c in range(5) for d in range(5) for e in range(2) if a+b+c+d <= 4]
	mat_list=[]
	for tuple in L:
		Y=M.subs(r==tuple[0],w==tuple[1],x==tuple[2],y==tuple[3],u==tuple[4])
		if spec_rad(Y)==6-s:
			if is_strongly_connected(Y):
				mat_list.append(Y)
	return(mat_list)

def tre(s):
	P=Matrix([[0,0,1,0],[1,0,0,0],[0,1,0,0],[0,0,0,1]])
	fin_list=[]
	mlist=uv_evals(s)
	for mat in mlist:
		fm=mpoly(mat,P,s)
		if fm in lam_polys(fm.coefficient(t^(4*s-1)),P,s):
			fin_list.append(mat)
	return(clear_dupes(fin_list))

###############################################################
#Two two-cycle
#Requires mathematica
###############################################################


def tt_gamcheck(s):
	P=Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
	L=[(i,j) for i in range(2*s+1) for j in range(i+1) if i+j<=2*s]
	for tup in L:
		lam=2*(tup[0]+tup[1])
		bet=tup[0]^2+tup[1]^2+4*tup[0]*tup[1]
		polys=lam_polys(-lam,P,s)
		gam_list=[bet-poly.coefficient(t^2) for poly in polys]
		max_gam=max(gam_list)
		print(tup,max_gam)
		#for poly in polys:
			#cp=poly.coefficient(t^2)
			#print(tup,lam,bet,bet-cp)

#We define the invariant z as z = (b-1)^2 + (h-1)^2 + 2*(c*e+d+f)
#Note that z is just gamma+2

#In the following routines, a==u, g==v
#tuples are (b,h,c,d,e,f)
#when ce=df=0 we have two cases in order to be strongly connected: either c=f=0 or d=e=0

def get_vars(tuple):
	all_vars=[x,y]
	tvars = [v for v in all_vars if v in tuple]
	return(tvars)

def get_divisors(n):
	div=divisors(n)
	L=[(d,n/d) for d in div]
	return(L)

def tbuilder(z):
	var('x y')
	L=[(b,h,u,v) for b in range(z+2) for h in range(z+2) for u in range(z+1) for v in range(z+1) if (b-1)^2+(h-1)^2+2*(u+v)==z]
	newL=[]
	for S in L:
		bvar,hvar,xvar,yvar=S
		#here x=c*e
		if xvar==0:
			xpairs=[[0,x],[x,0]]
		else:
			xpairs=get_divisors(xvar)
		#here y=d*f
		if yvar==0:
			ypairs=[[0,y],[y,0]]
		else:
			ypairs=get_divisors(yvar)
		newS = [(bvar,hvar,X[0],Y[0],X[1],Y[1]) for X in xpairs for Y in ypairs]
		newL = newL + newS
	return(newL)

def reducer(u,v,polys,z,s):
	var('a b c d e f g h x y')
	M=Matrix([[a,b,c,d],[b,a,d,c],[e,f,g,h],[f,e,h,g]])                       
	P=Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
	mp=mpoly(M,P,s)    
	tuples=tbuilder(z)
	all_solns=[]
	for L in tuples:
		p=mp.subs(a==u,g==v,b==L[0],h==L[1],c==L[2],d==L[3],e==L[4],f==L[5])
		psub=p.subs(t==1)
		tvars = get_vars(L)
		if len(tvars)>0:
			print(L)
			print(psub.factor())
			subsolve=psub.solve(tvars[0])
			print(subsolve)
			#for cond in subsolve:
				#print(p.subs(cond))

def zsolver(u,v,polys,z,s):
	m=mathematica
	var('a b c d e f g h x y')
	M=Matrix([[a,b,c,d],[b,a,d,c],[e,f,g,h],[f,e,h,g]])                       
	P=Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
	mp=mpoly(M,P,s) 
	tuples=tbuilder(z)
	all_solns=[]
	for L in tuples:
		p=mp.subs(a==u,g==v,b==L[0],h==L[1],c==L[2],d==L[3],e==L[4],f==L[5])
		psub=p.subs(t==1)
		for q in polys:
			chk=p-q
			if chk.coefficient(t^2)==0:
				tvars = get_vars(L)
				if len(tvars)==2:
					ineqs=[i>0 for i in tvars]
					sys=m([chk.coefficients(t)[0][0]==0]+[psub==0]+ineqs)
					S=sys.Solve(tvars,'Integers')
					#sys=[chk.coefficients(t)[0][0]]+[psub]+ineqs. --Maple code
					#S=maple(sys).solve() -- Maple code
					if len(S)>0: 
						#print(S) 
						ml=m(L)
						subs=(ml.ReplaceAll(S)).sage()
						for sub in subs:
							all_solns.append(tuple(sub))
				if len(tvars)==1:
					ineqs=[i>=0 for i in tvars]
					sys=m([chk.coefficients(t)[0][0]==0]+[psub==0]+ineqs)
					S=sys.Solve(tvars,'Integers')
					#sys=[chk.coefficients(t)[0][0]]+[psub]+ineqs
					#S=maple('sys').solve() -- Maple code
					#print(p.subs(t==1)) -- Maple code
					if len(S)>0: 
						#print(S) 
						ml=m(L)
						subs=(ml.ReplaceAll(S)).sage()
						for sub in subs:
							all_solns.append(tuple(sub))
				if len(tvars)==0 and chk==0:
					all_solns.append(L)
				#if len(S)>0:
					#print(q)
					#ml=m(L)
					#subs=(ml.ReplaceAll(S)).sage()
					#for sub in subs:
						#all_solns.append(tuple(sub))
	return sorted(set(all_solns))

def check_solns(u,v,polys,solns,s):
	var('a b c d e f g h x y')
	M=Matrix([[a,b,c,d],[b,a,d,c],[e,f,g,h],[f,e,h,g]])                       
	P=Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
	mp=mpoly(M,P,s) 
	all_mats=[]
	for S in solns:
		p=mp.subs(a==u,g==v,b==S[0],h==S[1],c==S[2],d==S[3],e==S[4],f==S[5])
		if p in polys:
			N=M.subs(a==u,g==v,b==S[0],h==S[1],c==S[2],d==S[3],e==S[4],f==S[5])		
			all_mats.append(N)
	return sort_and_deduplicate(all_mats)

def twotwocycle(u,v,maxz,s):
	P=Matrix([[0,1,0,0],[1,0,0,0],[0,0,0,1],[0,0,1,0]])
	all_mats=[]
	polys=lam_polys(-2*(u+v),P,s)
	for z in range(maxz+1):
		solns=zsolver(u,v,polys,z,s)
		all_mats=all_mats+check_solns(u,v,polys,solns,s)
	return(clear_dupes(all_mats))

def total_list(s):
	if s==3:
		zlist=[(2,2,0),(2,1,1),(2,0,4),(1,1,6),(1,0,5),(0,0,8)]
	if s==4:
		zlist=[(3,1,0),(2,2,2),(2,1,3),(2,0,6),(1,1,8),(1,0,7),(0,0,10)]
	all_mats=[]
	for T in zlist:
		all_mats=all_mats+twotwocycle(T[0],T[1],T[2],s) 
	return(all_mats)