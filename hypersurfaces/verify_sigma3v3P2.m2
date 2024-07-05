restart

A = QQ[apply(compositions(3,3), i-> a_i)]

aa = vars(A)

load "invariant_aronhold.m2"
glF = aa**diff(aa, Aro);

-- the dimension of its kernel
10^2 - numcols(mingens ideal glF) + 1
