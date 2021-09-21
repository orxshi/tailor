//=INPUTS=====================================================

W=5;
H=5;
D=1;
mx=100;
my=100;
mz=1;
bump=1;
radius = H/2;

//============================================================

pLL = newp; Point(pLL) = {-W,-H,0}; // lower left
pUL = newp; Point(pUL) = {-W,H,0};  // upper left
pLR = newp; Point(pLR) = {W,-H,0};  // lower right
pUR = newp; Point(pUR) = {W,H,0};   // upper right

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

outerbc[0] = out[2];
outerbc[1] = out[3];
outerbc[2] = out[4];
outerbc[3] = out[5];

emptybc[0] = out[0];
emptybc[1] = sBase;

Physical Surface(2) = outerbc[];
Physical Surface(3) = emptybc[];
Physical Volume(4) = {1};

//============================================================

//cnt_x = -W+W/3;
//cnt_y = 0;
//cnt_z = 0;
//cnt = newp; Point(cnt) = {cnt_x,cnt_y,cnt_z};
//
//pcircle_N = newp; Point(pcircle_N) = {cnt_x,cnt_y+radius,0};
//pcircle_S = newp; Point(pcircle_S) = {cnt_x,cnt_y-radius,0};
//pcircle_E = newp; Point(pcircle_E) = {cnt_x+radius,cnt_y,0};
//pcircle_W = newp; Point(pcircle_W) = {cnt_x-radius,cnt_y,0};
//
//circle_NW = newl; Circle(circle_NW) = {pcircle_N, cnt, pcircle_W};
//circle_WS = newl; Circle(circle_WS) = {pcircle_W, cnt, pcircle_S};
//circle_SE = newl; Circle(circle_SE) = {pcircle_S, cnt, pcircle_E};
//circle_EN = newl; Circle(circle_EN) = {pcircle_E, cnt, pcircle_N};
//
//lcntN = newl; Line(lcntN) = {cnt, pcircle_N};
//lcntS = newl; Line(lcntS) = {cnt, pcircle_S};
//lcntW = newl; Line(lcntW) = {cnt, pcircle_W};
//lcntE = newl; Line(lcntE) = {cnt, pcircle_E};
//
//Line Loop(1) = {lcntN, circle_NW, -lcntW};
//Line Loop(2) = {lcntW, circle_WS, -lcntS};
//Line Loop(3) = {lcntS, circle_SE, -lcntE};
//Line Loop(4) = {lcntE, circle_EN, -lcntN};
//
//sNW = news; Plane Surface(sNW) = {1};
//sWS = news; Plane Surface(sWS) = {2};
//sSE = news; Plane Surface(sSE) = {3};
//sEN = news; Plane Surface(sEN) = {4};
//
//out[] = Extrude {0,0,D} {Surface{sNW,sWS,sSE,sEN}; Layers{mz}; Recombine;};
//
//Transfinite Line{circle_NW, circle_WS, circle_SE, circle_EN} = 10 Using Progression 1.0;
//Transfinite Line{lcntN} = 10 Using Progression 1.0;
//Transfinite Line{lcntS} = 10 Using Progression 1.0;
//Transfinite Line{lcntW} = 10 Using Progression 1.0;
//Transfinite Line{lcntE} = 10 Using Progression 1.0;
//
//Transfinite Surface{sNW,sWS,sSE,sEN}; Recombine Surface{sNW,sWS,sSE,sEN};
//
//Printf("out[0]: %g", out[0]);
//Printf("out[1]: %g", out[1]);
//Printf("out[2]: %g", out[2]);
//Printf("out[3]: %g", out[3]);
//Printf("out[4]: %g", out[4]);
//Printf("out[5]: %g", out[5]);
//Printf("out[6]: %g", out[6]);
//Printf("out[7]: %g", out[7]);
//Printf("out[8]: %g", out[8]);
//Printf("out[9]: %g", out[9]);
//Printf("out[10]: %g", out[10]);
//Printf("out[11]: %g", out[11]);
//Printf("out[12]: %g", out[12]);
//Printf("out[13]: %g", out[13]);
//Printf("out[14]: %g", out[14]);
//Printf("out[15]: %g", out[15]);
//Printf("out[16]: %g", out[16]);
//Printf("out[17]: %g", out[17]);
//Printf("out[18]: %g", out[18]);
//Printf("out[19]: %g", out[19]);
//
//interogbc[0] = out[3];
//interogbc[1] = out[8];
//interogbc[2] = out[13];
//interogbc[3] = out[18];
//
//emptybc[0] = sNW;
//emptybc[1] = sWS;
//emptybc[2] = sSE;
//emptybc[3] = sEN;
//emptybc[4] = out[0];
//emptybc[5] = out[5];
//emptybc[6] = out[10];
//emptybc[7] = out[15];
//
//Physical Surface(11) = interogbc[];
//Physical Surface(3) = emptybc[];
//Physical Volume(4) = {1,2,3,4};
