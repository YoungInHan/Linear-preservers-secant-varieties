restart

KK = QQ

T = KK[t_(0,0,0)..t_(2,2,2)]


ttt = vars(T)

load "invariant_strassen.m2"
-- Strassen's invariant is loaded;
-- it is a polynomial of degree 9 in 27 variables;
-- in the expansion in the monomial basis t_(i,j,k) it has 9216 terms;

lieF = ttt**diff(ttt,Str);

27^2 - numcols(mingens ideal lieF) + 1
-- the dimension of the preserver is 25
-- the expected dimension is 8+8+8+1 = dim G(C3,C3,C3)
