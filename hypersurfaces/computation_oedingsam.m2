restart


-- we set up a ring with variables vk_i playing the roles of basis vectors
-- and t_(i1,..,i5) playing the role of the coefficients of the tensor

V1 = QQ[v1_0..v1_1]
V2 = QQ[v2_0..v2_1]
V3 = QQ[v3_0..v3_1]
V4 = QQ[v4_0..v4_1]
V5 = QQ[v5_0..v5_1]

-- weight for the SL2 action:
wts = toList apply((0,0,0,0,0)..(1,1,1,1,1), iii-> (
	{1}| apply(toList iii, j-> (-1)^j)))

T = QQ[t_(0,0,0,0,0)..t_(1,1,1,1,1), Degrees =>  wts]


TABC = T**V1**V2**V3**V4**V5

ttt = sub(vars(T),TABC)

vv1 = sub(vars(V1),TABC)
vv2 = sub(vars(V2),TABC)
vv3 = sub(vars(V3),TABC)
vv4 = sub(vars(V4),TABC)
vv5 = sub(vars(V5),TABC)


-- the zero weight space for the action of SL x SL x SL x SL x SL
ttt6wt0 = sub(basis({6,0,0,0,0,0},T),TABC);

-- the degree 6 invariant is skew with respect to the action of S5
-- which permutes the five factors. Skew-symmetrization reduces the 
-- dimension of the search space
permOpts = ppp -> toList apply( (0,0,0,0,0)..(1,1,1,1,1), iii->( t_iii => t_(toSequence(iii_ppp))))
permSgn = ppp -> (det( (id_(QQ^(#ppp)))_ppp))

U = sum(permutations(5), ppp -> (
	print ppp;
	permSgn(ppp) * sub(ttt6wt0, permOpts(ppp))));
U = mingens ideal U;


-- we interpolate on the resulting 25 generators
evalOpts = Tt -> (
    toList apply((0,0,0,0,0)..(1,1,1,1,1), iii-> (
	    t_iii => diff( v1_(iii_0)*v2_(iii_1)*v3_(iii_2)*v4_(iii_3)*v5_(iii_4),Tt))))


A = matrix apply(25, j-> (
	print j;
	Tt = sum(5, j-> (
		sub(random(1,V1),TABC) * 
		sub(random(1,V2),TABC) * 
		sub(random(1,V3),TABC) * 
		sub(random(1,V4),TABC) * 
		sub(random(1,V5),TABC) ));
	flatten entries sub(U,evalOpts(Tt))));

-- compute the kernel to obtain the desired invariant
F = (U * gens ker A)_(0,0);

"invariant_oedingsam.m2" << "OSinv = ";
"invariant_oedingsam.m2" << toString(F);
"invariant_oedingsam.m2" << ";" << close


