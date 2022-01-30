#############################################################################
##  Shamelessly ripped and hacked from qpa source files
##
##
##  The primary function herein is QuiverHilbert(M,d,spdeg)
##  which inputs an adjacency matrix M and computes the d coefficient Matrix-valued
##  Hilbert series corresponding to a random (homogeneous) superpotential
##  of degree spdeg.
LoadPackage("qpa");

# create random (rational) coefficients
RandRat := function( n )

	local b,rat;

	b:=0; 
	while (b=0) do 
		b:=Random([-n..n]); 
	od;
		

	rat := Random([-n..n]/b); 
	return(rat);
end;

UnWalk:=function( w ) 
    local path,len,i;

    path := w[1];
    len := Length(w);

    for i in [2..len] do
	path := path*w[i];
    od;

    return path;
end;

## Generalization of PathsOfLengthTwo
## Inputs a quiver and an integer (length)
## Outputs a list of paths of the given length
PathsOfLength:=
    function( q, len ) 
 
    local paths, new, p, t, out, a;

    # Recursively get paths of length len-1
    if (len = 1) then
	new := ArrowsOfQuiver(q);
    elif (len = 2) then
	new := PathsOfLengthTwo(q);
    else
	paths := PathsOfLength(q,len-1);
    
	new := [];    
  
    	for p in paths do
		# t holds the target of p
		t := TargetOfPath(p);
        	# out is all the arrows out of t
        	out := OutgoingArrowsOfVertex(t); 
        	for a in out do
			Add(new,p*a);
		od;
    	od;
    fi;

    return new;
end; 

CyclicPaths:=function( q, spdeg ) 
 
    local cycles,paths,np,p,s,t,g,perms,w,i,obj;

    # c will hold the list of all cycles of length len
    cycles := [];
  
    # paths is a list of all paths of length len
    paths := PathsOfLength(q,spdeg);
    np := Length(paths);

    while (np > 0) do  
	p := paths[1];
	Remove(paths,1);

	# s holds the source of p and t holds the target
	s := SourceOfPath(p);
	t := TargetOfPath(p);
	# first check that they are the same and if so add to cycles
	# then remove all permutations from paths
	if (s=t) then
		Add(cycles,p);

		# make a list of all cyclic permutations of p
		g := GeneratorsOfGroup(SymmetricGroup(spdeg))[1];
		perms := [];
		w := WalkOfPath(p);
		for i in [1..spdeg-1] do
			Add(perms,UnWalk(Permuted(w,g^i)));
		od;

		# find each element of perms and remove from paths
		for obj in perms do 
			if (obj in paths) then
				Remove(paths,Position(paths,obj));	
			fi;		
		od;
	fi;

	np := Length(paths);
    od;

    Sort(cycles);
    return cycles;
end;

##  This function return a list of all paths 
##  determined by derivations of 3 cycles
##  by all paths
ZeroPaths:=
    function( q ) 
 
    local pa, zero, gens, numg, numv, rels, i,j;

    # zero element and generators
    pa := PathAlgebra(Rationals,q);
    zero := Zero(pa);
    gens := GeneratorsOfAlgebra(pa);
    numg := Length(gens);

    # Get number of trivial paths
    numv := NumberOfVertices(q);	  

    # Calculate zero relations:
    rels := [];
    for i in [numv+1..numg] do
      for j in [numv+1..numg] do
        # Check if relation is zero
	if (gens[i]*gens[j] = zero) then
		# Add relation and adjust indexing to account for trivial paths
          	Add(rels,[[[i-numv,j-numv]],[1]]);
	fi;
      od;
    od;

    return rels;
end;

# This function takes a quiver and a degree
# and computes a random superpotential of that degree
# and returns the relations determined by taking
# n left derivations.
Superpotential:=function( q, spdeg, n ) 
 
    local arrows,narr,nrel,cyc,tot,cof,sup,g,C,i,j,w,wlist,k,l,per,lts,pos;

    # underlying quiver
    arrows := ArrowsOfQuiver(q);
    narr := NumberOfArrows(q);
 
    # three cycles to derivate
    cyc := CyclicPaths(q,spdeg);
    tot := Length(cyc);

    # create random (rational) coefficients for the cyclic terms
    cof := [];
    for i in [1..tot] do
	Add(cof,RandRat(50));
    od;

    # Get relations from superpotential
    sup := []; 
    g := GeneratorsOfGroup(SymmetricGroup(spdeg))[1];

    # Get possible derivations tuples and create
    # an empty list to store relations
    C := UnorderedTuples([1..narr],n);	
    nrel := Length(C);
    for i in [1..nrel] do sup[i] := [[],[]]; od;

    # for each cyclic path, first replace path notation with
    # list notation then permute and take derivations
    for j in [1..tot] do
   	w:=WalkOfPath(cyc[j]);
	wlist := [];
	# Replace path notation with list notation
	for k in [1..spdeg] do 
		for l in [1..narr] do
			if arrows[l]=w[k] then
				wlist[k] := l;
				break;
			fi;
		od;
	od;    	
    	
	# permute and then chop off first n elements and add to sup
    	for i in [0..spdeg-1] do 
		# permute cyclic term
		per := Permuted(wlist,g^i);
		
		# get leading terms, save in lts, then chop off
		lts := [];
		for k in [1..n] do
			Add(lts,per[1]);
			Remove(per,1);
		od;

		# find lts in C
		Sort(lts);
		pos := Position(C, lts);
		
		# add remaining relation to sup and its corresponding coefficient
		Add(sup[pos][1],per); Add(sup[pos][2],cof[j]); 
	od;	
    od;
    return sup;
end;

## Puts all the quiver relations in GBNP format
## then computes a basis of degree d
QuiverTerms:=function( M, d, spdeg ) 
 
    local q,sup,rels,srel,fgens,wt,i,A,GB,dims,BT,der;

    #Set number of derivations
    if spdeg=2 then 
	der := 0;
    else
	der := 1;
    fi;		

    #Get quiver relations
    q := Quiver(M);
    sup := Superpotential(q,spdeg,der);
    rels := ZeroPaths(q);
    for srel in sup do
	Add(rels,srel);
    od;
	
    # generators of path algebra ignoring trivial paths
    fgens := Length(ArrowsOfQuiver(q));

    #Weights of generators - all 1 at this time
    wt:=[];
    for i in [1..fgens] do
	Add(wt,1);
    od;

    #The free algebra and its zero paths
    A:=FreeAssociativeAlgebraWithOne(Rationals,fgens);

    #Call GBNP
    GBNP.ConfigPrint(A);
    GB:=SGrobnerTrunc(rels, d, wt);
    dims:=DimsQATrunc(GB, d, wt);
    Print(dims);
    Print("\n");
    BT:=BaseQATrunc(GB, d, wt);
    
    #for degpart in BT do for mon in degpart do PrintNP([[mon],[1]]); od; od;

    return(BT[d+1]);
end;

# Primary function.
# Inputs an adjacency matrix M,
# degree d for term in Hilbert series
# and degree spdeg of the corresponding superpotential.
QuiverHilbert:=function( M, d, spdeg ) 
 
    local rdim,cdim,mons,T,tot,i,j,m,k,L,N,term,first,last,out,u,v;
   
    if (not IsMatrix(M)) then
        Print("QuiverHilbert(): argument is not a matrix\n");
        return fail;
    fi;

    rdim := DimensionsMat(M)[1];
    cdim := DimensionsMat(M)[2];

    if (rdim <> cdim) then
	Print("QuiverHilbert(): argument is not a square matrix\n");
        return fail;
    fi;

    #Get number of monomials of degree d
    mons := QuiverTerms(M,d,spdeg);
    
    #Create matrix of lists with arrows corresponding to entries
    T := 0*IdentityMat(rdim);
    tot := 1;
    for i in [1..rdim] do
   	for j in [1..rdim] do	    
		m := M[i][j];
		T[i][j] := [];
		for k in [tot..tot+m-1] do
			Add(T[i][j],k);
			tot := tot + 1;
		od;
	od;
     od;

    #Matrix to count terms
    L := [];
    N := 0*IdentityMat(rdim);

    #Find paths coming out of each vertex
    for term in mons do
        first := term[1]; 
	last := term[d];
	
	out := false;

	for i in [1..rdim] do
	for j in [1..rdim] do
		if (first in T[i][j]) then 
			for u in [1..rdim] do
			for v in [1..rdim] do
				if (last in T[u][v]) then N[i][v] := N[i][v]+1; fi;
			od;od;
		fi;
	od;od;
    od;

    return N;
end;

# To determine superpotential relations
# with symbolic coefficients

Symbsp:=function( q, spdeg, n, coeffs ) 
 
    local arrows,narr,nrel,cyc,tot,sup,g,C,i,j,w,wlist,k,l,per,lts,pos;

    # underlying quiver
    arrows := ArrowsOfQuiver(q);
    narr := NumberOfArrows(q);
 
    # three cycles to derivate
    cyc := CyclicPaths(q,spdeg);
    tot := Length(cyc);

    # create random (rational) coefficients for the cyclic terms
    #coeffs := ["c1","c2","c3","c4","c5","c6","c7","c8","c9"];

    # Get relations from superpotential
    sup := []; 
    g := GeneratorsOfGroup(SymmetricGroup(spdeg))[1];

    # Get possible derivations tuples and create
    # an empty list to store relations
    C := UnorderedTuples([1..narr],n);	
    nrel := Length(C);
    for i in [1..nrel] do sup[i] := [[],[]]; od;

    # for each cyclic path, first replace path notation with
    # list notation then permute and take derivations
    for j in [1..tot] do
   	w:=WalkOfPath(cyc[j]);
	wlist := [];
	# Replace path notation with list notation
	for k in [1..spdeg] do 
		for l in [1..narr] do
			if arrows[l]=w[k] then
				wlist[k] := l;
				break;
			fi;
		od;
	od;    	
    	
	# permute and then chop off first n elements and add to sup
    	for i in [0..spdeg-1] do 
		# permute cyclic term
		per := Permuted(wlist,g^i);
		
		# get leading terms, save in lts, then chop off
		lts := [];
		for k in [1..n] do
			Add(lts,per[1]);
			Remove(per,1);
		od;

		# find lts in C
		Sort(lts);
		pos := Position(C, lts);
		
		# add remaining relation to sup and its corresponding coefficient
		Add(sup[pos][1],per); Add(sup[pos][2],coeffs[j]); 
	od;	
    od;
    return sup;
end;

#Returns the paths of a length n in a quiver q
#between the given vertices i,j
PathsBtw:=function( q, n, i, j)

	local P,L,p,V;

	V:=VerticesOfQuiver(q);

	P:=PathsOfLength(q,n);
	L:=[];
	for p in P do
		if (SourceOfPath(p)=V[i] and TargetOfPath(p)=V[j]) then 
			Add(L,p); 
		fi;
	od;
	return L;
end;

#Given a walk of a path w in a 
#returns all subpaths of length min m and max M
SubPaths:=function( w, m, M)

	local len,subs,k,j;

	len := Number(w);

	subs:=[];
	for k in [m..M] do
		for j in [1..len-k+1] do
			Add(subs,w{[j..j+k-1]});
		od;
	od;
	return(subs);
end;

#Returns the paths of length n in a quiver q
#modulo the relations R between vertices i,j
PathsMod:=function( q, n, R, i, j ) 
 
    	local arrows,verts,P,paths,p,w,subs;

    	# underlying quiver
    	arrows := ArrowsOfQuiver(q);
    	verts := VerticesOfQuiver(q);
   	P:=PathsBtw(q,n,i,j); 

	paths := [];

    	for p in P do
		w:=WalkOfPath(p);
		subs:=SubPaths(w,1,n);

		if (Intersection(subs,R)=[]) then
			Add(paths,p);
		fi;
	od;

	return paths;
end;

#Compiles PathsMod into a matrix
FullPathsMod:=function(q,n,R)

	local i,j,v,T;

	v := Number(VerticesOfQuiver(q));
	T := 0*IdentityMat(v);

	for i in [1..v] do
		for j in [1..v] do
			T[i][j] := Number(PathsMod(q,n,R,i,j));
	od; od;

	return T;
end;

#Given a group G and the character of a representation V,
#this function computes the McKay quiver corresponding to V.
#The result is given as the adjacency quiver of the McKay quiver.
mckayQuiver :=function(G,chV) 

	local t,m,n,x,i,j,M,T;

	if (not IsGroup(G)) then
        	Print("mckayQuiver(): argument is not a group\n");
        	return fail;
	fi;

	t := CharacterTable(G);
	n := Length(ConjugacyClasses(G));
	
	#character table as a matrix
	m := List(Irr(t),ValuesOfClassFunction);
	T := TransposedMat(m);

	x := 0*IdentityMat(n);
	#Compute the character of the tensor product
	#of V and V_i for each irreducible V_i
	for i in [1..n] do
	for j in [1..n] do
		x[i][j] := chV[j]*m[i][j];
	od;od;

	M := (Inverse(T)*TransposedMat(x));

	return M;
end;

#The function prints those small groups that have c conjugacy classes.
#Works for c at most 8.
cclass := function(c) 

	local L,n,j,A,nA,i,G,l,V;

	L:=[1,2,6,42,1806,3263442,10650056950806,113423713055421844361000442,
12864938683278671740537145998360961546653259485195806];
	n:=L[c];

	V:=[];

    for j in [1..n] do
	A := AllSmallGroups(j);
	nA := Length(A);

	for i in [1..nA] do
		G:=A[i];
		l := Length(ConjugacyClasses(G));
		if (l=c) then Add(V,[j,i]); fi;
	od;
    od;

return(V);
end;

#Prints the mckay quivers for groups G with character x
#of dimension d
mq := function(G,d)

	local I,nI,i,dI,A,nA,a,x,k;

	I := Irr(G);
	nI := Length(I);

	dI:=[];
	for i in [1..nI] do 
		Add(dI,I[i][1]);
	od;

	for i in [1..d] do
		A := UnorderedTuples([1..nI],i); 
		nA:=Length(A);
		
		for a in A do
			x:=Sum(a,k->I[k]);
			if x[1]=d then
				Print("x=", a, ", ", "M=", mckayQuiver(G,x), "\n");
			fi;
		od;
	od;
	
end;


