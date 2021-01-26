//=INPUTS=====================================================

H=7;
D=1;
mx=20;
my=20;
mz=1;
bump=1;

//============================================================

pLL = newp; Point(pLL) = {-H,-H,0}; // lower left
pUL = newp; Point(pUL) = {-H,H,0};  // upper left
pLR = newp; Point(pLR) = {H,-H,0};  // lower right
pUR = newp; Point(pUR) = {H,H,0};   // upper right

lL = newl; Line(lL) = {pLL,pUL};
lB = newl; Line(lB) = {pLR,pLL};
lR = newl; Line(lR) = {pUR,pLR};
lT = newl; Line(lT) = {pUL,pUR};

Transfinite Line{lL} = my Using Bump bump;
Transfinite Line{lR} = my Using Bump bump;
Transfinite Line{lB} = mx Using Bump bump;
Transfinite Line{lT} = mx Using Bump bump;

Line Loop(1) = {lL, lT, lR, lB};
sBase = news; Plane Surface(sBase) = {1};

Transfinite Surface{sBase}; Recombine Surface{sBase};

