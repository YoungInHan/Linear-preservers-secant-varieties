restart

-- the calculation to generate the invariant can be performed over a finite field
-- the resulting invariant has coefficient +/-1 and it can be interpreted over the integers
-- one can then verify that it is the correct invariant by explicitly evaluating on 
-- low rank tensors

KK = ZZ/101
KK = QQ

-- we set up a ring with variables a_i,b_j,c_k playing the roles of basis vectors
-- and t_(i,j,k) playing the role of the coefficients of the tensor
A = KK[a_0..a_2]
B = KK[b_0..b_2]
C = KK[c_0..c_2]

T = KK[t_(0,0,0)..t_(2,2,2), Degrees =>  toList({1,0,0,0}..{1,2,2,2})]

TABC = T**A**B**C

ttt = sub(vars(T),TABC)

aa = sub(vars(A),TABC)
bb = sub(vars(B),TABC)
cc = sub(vars(C),TABC)

-- the zero weight space for the action of SL x SL x SL 
ttt9wt0 = sub(basis({9,9,9,9},T),TABC);

-- Strassen's invariant is skew with respect to the action of S3
-- which permutes the three factors. Skew-symmetrization reduces the 
-- dimension of the search space
permOpts = ppp -> toList apply( (0,0,0)..(2,2,2), iii->( t_iii => t_(toSequence(iii_ppp))))
permSgn = ppp -> (det( (id_(ZZ^(#ppp)))_ppp))

U = sum(permutations(3), ppp -> (
	print ppp;
	permSgn(ppp) * sub(ttt9wt0, permOpts(ppp))));
U = mingens ideal U;

-- invariancy for the action of SL x SL x SL implies invariancy
-- under the action of the Weyl group S3 x S3 x S3.
permWeyl1 = ppp -> toList apply( (0,0,0)..(2,2,2), iii->( t_iii => t_(toSequence(ppp_(iii_0),iii_1,iii_2))))
permWeyl2 = ppp -> toList apply( (0,0,0)..(2,2,2), iii->( t_iii => t_(toSequence(iii_0,ppp_(iii_1),iii_2))))
permWeyl3 = ppp -> toList apply( (0,0,0)..(2,2,2), iii->( t_iii => t_(toSequence(iii_0,iii_1,ppp_(iii_2)))))

U1 = sum(permutations(3), ppp -> (
	print ppp;
	permSgn(ppp) * sub(U, permWeyl1(ppp))));
U1 = mingens ideal U1;
U2 = sum(permutations(3), ppp -> (
	print ppp;
	permSgn(ppp) * 	sub(U1, permWeyl2(ppp))));
U2 = mingens ideal U2;
U3 = sum(permutations(3), ppp -> (
	print ppp;
	permSgn(ppp) * 	sub(U2, permWeyl3(ppp))));
U3 = mingens ideal U3;



-- we interpolate on the resulting 1225 generators
evalOpts = Tt -> (
    toList apply((0,0,0)..(2,2,2), iii-> (
	    t_iii => diff( a_(iii_0)*b_(iii_1)*c_(iii_2),Tt))))


A = matrix apply((numcols(U3), j-> (
	print j;
	Tt = sum(4, k-> (
		sub(random(1,A),TABC) * 
		sub(random(1,B),TABC) * 
		sub(random(1,C),TABC)  ));
	flatten entries sub(U3,evalOpts(Tt)))));

-- compute the kernel to obtain the desired invariant
A = sub(A,KK);
kA = gens ker A;

kA -- observe that only the top 17 entries of kA are nonzero
F = (U3_{0..16} * kA^{0..16})_(0,0);

"strassen_invariant.m2" << "Str = ";
"strassen_invariant.m2" << toString(F);
"strassen_invariant.m2" << ";" << close

load "strassen_invariant.m2"

