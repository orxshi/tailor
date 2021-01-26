//=INPUTS=====================================================

H=20;
I=H/2;
D=1;
mx=10;
my=10;
mz=1;
bump=1;
originx=-50;
originy=0;

//============================================================

oLL = newp; Point(oLL) = {originx-H,originy-H,0}; // lower left
oUL = newp; Point(oUL) = {originx-H,originy+H,0};  // upper left
oLR = newp; Point(oLR) = {originx+H,originy-H,0};  // lower right
oUR = newp; Point(oUR) = {originx+H,originy+H,0};   // upper right

iLL = newp; Point(iLL) = {originx-I,originy-I,0}; // lower left
iUL = newp; Point(iUL) = {originx-I,originy+I,0};  // upper left
iLR = newp; Point(iLR) = {originx+I,originy-I,0};  // lower right
iUR = newp; Point(iUR) = {originx+I,originy+I,0};   // upper right

iUL_oUL = newl; Line(iUL_oUL) = {iUL,oUL};
oUL_oUR = newl; Line(oUL_oUR) = {oUL,oUR};
oUR_iUR = newl; Line(oUR_iUR) = {oUR,iUR};
iUR_iUL = newl; Line(iUR_iUL) = {iUR,iUL};

oLL_oUL = newl; Line(oLL_oUL) = {oLL,oUL};
iUL_iLL = newl; Line(iUL_iLL) = {iUL,iLL};
iLL_oLL = newl; Line(iLL_oLL) = {iLL,oLL};

iLL_iLR = newl; Line(iLL_iLR) = {iLL,iLR};
iLR_oLR = newl; Line(iLR_oLR) = {iLR,oLR};
oLR_oLL = newl; Line(oLR_oLL) = {oLR,oLL};

oUR_oLR = newl; Line(oUR_oLR) = {oUR,oLR};
iLR_iUR = newl; Line(iLR_iUR) = {iLR,iUR};

Transfinite Line{oUL_oUR} = my;
Transfinite Line{oUR_oLR} = my;
Transfinite Line{oLR_oLL} = my;
Transfinite Line{oLL_oUL} = my;

Transfinite Line{iUR_iUL} = my;
Transfinite Line{iLR_iUR} = my;
Transfinite Line{iLL_iLR} = my;
Transfinite Line{iUL_iLL} = my;

Transfinite Line{iUL_oUL} = my Using Progression bump;
Transfinite Line{iLR_oLR} = my Using Progression bump;
Transfinite Line{-oUR_iUR} = my Using Progression bump;
Transfinite Line{iLL_oLL} = my Using Progression bump;

Line Loop(1) = {iUR_iUL, iUL_oUL, oUL_oUR, oUR_iUR};
Line Loop(2) = {-iUL_oUL, iUL_iLL, iLL_oLL, oLL_oUL};
Line Loop(3) = {-iLL_oLL, iLL_iLR, iLR_oLR, oLR_oLL};
Line Loop(4) = {-iLR_oLR, iLR_iUR, -oUR_iUR, oUR_oLR};

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};

Transfinite Surface{1}; Recombine Surface{1};
Transfinite Surface{2}; Recombine Surface{2};
Transfinite Surface{3}; Recombine Surface{3};
Transfinite Surface{4}; Recombine Surface{4};
