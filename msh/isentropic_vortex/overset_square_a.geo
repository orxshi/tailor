//=INPUTS=====================================================

Wbg = 20;
Wdiff = 6;
W=1;
H=3;
D=1;
mx=50;
my=150;
mz=1;
bump=1;
radius = H/2;

//============================================================

pLL = newp; Point(pLL) = {-W-Wbg+W+Wdiff,-H,0}; // lower left
pUL = newp; Point(pUL) = {-W-Wbg+W+Wdiff,H,0};  // upper left
pLR = newp; Point(pLR) = {W-Wbg+W+Wdiff,-H,0};  // lower right
pUR = newp; Point(pUR) = {W-Wbg+W+Wdiff,H,0};   // upper right

lL = newl; Line(lL) = {pLL,pUL};
lB = newl; Line(lB) = {pLR,pLL};
lR = newl; Line(lR) = {pUR,pLR};
lT = newl; Line(lT) = {pUL,pUR};

Line Loop(1) = {lL, lT, lR, lB};
sBase = news; Plane Surface(sBase) = {1};

Transfinite Line{lL,lR} = my Using Progression 1.0;
Transfinite Line{lB,lT} = mx Using Progression 1.0;
Transfinite Surface{sBase}; Recombine Surface{sBase};

out[] = Extrude {0,0,D} {Surface{sBase}; Layers{mz}; Recombine;};

Printf("out[0]: %g", out[0]);
Printf("out[1]: %g", out[1]);
Printf("out[2]: %g", out[2]);
Printf("out[3]: %g", out[3]);
Printf("out[4]: %g", out[4]);
Printf("out[5]: %g", out[5]);

interogbc[0] = out[2];
interogbc[1] = out[3];
interogbc[2] = out[4];
interogbc[3] = out[5];

emptybc[0] = out[0];
emptybc[1] = sBase;

Physical Surface(11) = interogbc[];
Physical Surface(3) = emptybc[];
Physical Volume(4) = {1};
