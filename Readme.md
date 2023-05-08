Here is a description of the relevant files for our Sage computations. See below for the (mostly) deprecated routines in Maple.

Files: cypack.sage

Descriptions.

cypack.sage This is the primary set of functions.

spec_rad(M) - returns the spectral radius of a matrix M.

is_strongly_connected(M) - checks whether a matrix M is strongly connected. 

is_normal(M) - checks whether a matrix M is normal.

mpoly(M,P,s) - returns the determinant of the matrix polynomial for M with permutation matrix P and superpotential of degree s.

hilb_series(M,P,n,s) - returns the Hilbert series for a CY matrix, up to degree n, with data as in mpoly(M,P,s).

sort_cycs(deg) - returns cyclotomic polynomials of degree deg such that (1-t)^3 is a factor.

clear_dupes(mlist) - cleans a list of matrices removing those that are equivalent up to permuting the vertices.

mutate(M,v) - computes the adjacency matrix of a the quiver with adjacency matrix M after mutation at vertex v.

CYpoly(M,d) - returns the matrix polynomial associated to a quiver supporting a CY3 algebra with a superpotential of degree d.


___________________________________________________________________
The following is the old Readme for Maple files. Most of these are now incorporated in the cypack.sage file. The key exception is the mfind routine which will replaced in future work. 
Here is a description of the relevant files for our Maple computations.
Files: cypack.mw / cypack.m ListSpdeg3.mw ListSpdeg4.mw
Descriptions.
cypack.mw / cypack.m This is the primary set of functions.
ListSpdeg3.mw / ListSpdeg4.mw Runs the search procedure for degree 3/4 superpotentials and produces all "good" matrices. Many of these are not actually good but Maple cannot judge them.
IsStronglyConnected(M) - checks whether a matrix M is strongly connected. This is implemented in mforms.mw to generate generic forms of potential matrices.
IsNormal(M) - checks whether a matrix M is normal.
CYpoly(M,d) - returns the matrix polynomial associated to a quiver supporting a CY3 algebra with a superpotential of degree d.
mpoly(M,d) - returns the determinant of CYpoly(M,d).
CY(M,d) - computes the roots of mpoly(M,d).
hseries(M,n,d) - returns a list of the first n terms of the matrix-valued Hilbert series associated to a quiver supporting a CY3 algebra with a superpotential of degree d.
thilb(M,n,d) - returns a list of the first n terms in hseries(M,n,d) after summing all entries in each term.
ph(M,n,d) - prints hseries(M,n,d). Helpful if you want the output as opposed to the list.
hpoly(M,n,d) - returns hseries(M,n,d) written as a polynomial in t up to degree n.
IsGood(M,d) - checks whether all roots of mpoly(M,d) lie on the unit circle.
first(M,d) - evaluates the derivative of mpoly(M,d) at t=1.
GK(M,d) - computes the GK dimension of a vacualgebra A living on the quiver with adjacency matrix M with superpotential of degree d, assuming such a thing exists.
mutate(M,v) - computes the adjacency matrix of a the quiver with adjacency matrix M after mutation at vertex v.
the next set of functions are used to construct generic forms for matrices in our classification.
mfind(u,v,w,gam,d) - first grabs all generic matrix forms with diagonal entries u,v,w and gamma=gam. For each possible form, the matrix polynomial is computed and compared to potentially valid products of cyclotomic polynomials. If a solution is found for indeterminates in generic forms then this solution is returned.
