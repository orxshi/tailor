//=INPUTS=====================================================

H=20;
//mx=20;
//my=20;
//mz=1;
bump=1;
depth = H*2; // extrusion length

//============================================================

pLL = newp; Point(pLL) = {-H,-H,-H}; // lower left
pUL = newp; Point(pUL) = {-H,H,-H};  // upper left
pLR = newp; Point(pLR) = {H,-H,-H};  // lower right
pUR = newp; Point(pUR) = {H,H,-H};   // upper right

lL = newl; Line(lL) = {pLL,pUL};
lB = newl; Line(lB) = {pLR,pLL};
lR = newl; Line(lR) = {pUR,pLR};
lT = newl; Line(lT) = {pUL,pUR};

//Transfinite Line{lL} = my Using Bump bump;
//Transfinite Line{lR} = my Using Bump bump;
//Transfinite Line{lB} = mx Using Bump bump;
//Transfinite Line{lT} = mx Using Bump bump;

Line Loop(1) = {lL, lT, lR, lB};
sBase = news; Plane Surface(sBase) = {1};

//Transfinite Surface{sBase};
//Recombine Surface{sBase};

//out[] = Extrude {0,0,depth} {Surface{sBase}; Layers{1}; Recombine;};
out[] = Extrude {0,0,depth} {Surface{sBase};};

center = newp;
Point(center) = {0,0,0};

lc = 10;
Field[1] = Distance;
Field[1].NodesList = {center};
Field[2] = MathEval;
Field[2].F = Sprintf("F1/1 + %g", lc / 1000);
Background Field = 2;

Physical Surface(9) = {out[2], out[3], out[4], out[5], sBase, out[0]};
Physical Volume(4) = {out[1]};
