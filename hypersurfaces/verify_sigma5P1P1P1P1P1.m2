restart

T = QQ[t_(0,0,0,0,0)..t_(1,1,1,1,1), Degrees =>  toList({1,0,0,0,0,0}..{1,1,1,1,1,1})]


ttt = vars(T)

load "invariant_oedingsam.m2"

lieF = ttt**diff(ttt,OSinv);

32^2 - numcols(mingens ideal lieF) + 1
