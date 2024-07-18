restart

-- we set up a ring with variables x_i playing the role of basis vectors
-- and a_(i,j,k) playing the role of the coefficients of the homogeneous polynomial
R = QQ[x_0..x_2]
A = QQ[apply(compositions(3,3), i-> a_i)]

AR = A**R

xx = sub(vars(R),AR)
aa = sub(vars(A),AR)

fgen = 6*sum(compositions(3,3), cc -> (
	a_cc*product(3, j-> (1/(cc_j)! * x_j^(cc_j)))))

-- the koszul map
kos = (res(ideal(xx))).dd_2 * matrix{{0,0,1},{0,-1,0},{1,0,0}}

semiFlat = diff(kos,diff(xx,diff(transpose xx, fgen)))
-- semiFlat is skew-symmetric, all its 8x8 pfaffians coincide and 
-- Aronhold invariant is any of such pfaffians
F = 1/1296* (pfaffians(8, semiFlat))_0


"invariant_aronhold.m2" << "Aro = ";
"invariant_aronhold.m2" << toString(F);
"invariant_aronhold.m2" << ";" << close


