c=1; // chord length
t=12 * c / 100; // thickness -- NACA0012
nX=1000; // number of 'x' points
depth=1000; // extrusion length
transX=30; //100
transY=30; //100
transWake=30; //100
stretchWake=0.95;
major=30;
minor=30;
originx=0;
originy=0;
nz=1;
stretchVert=0.9;

cenPx = originx;
cenPy = originy;

xInc = c / nX;
C1 = 0.2969;
C2 = 0.1260;
C3 = 0.3516;
C4 = 0.2843;
C5 = 0.1036;

x[0] = 0;
For i In {1:nX}
	x[i] = i * xInc;
EndFor

y[0]  = 0;
For i In {1:nX-1}
	y[i] = t * c / 0.2 * ( C1 * (x[i]/c)^(1/2) - 
			       C2 * (x[i]/c)^(1)   -
			       C3 * (x[i]/c)^(2)   +
			       C4 * (x[i]/c)^(3)   -
			       C5 * (x[i]/c)^(4)   );
EndFor
y[nX] = 0;

For i In {0:nX}
	pFoilUpper[i] = newp;
	Point( pFoilUpper[i] ) = {x[i],y[i],0};
EndFor

Translate {cenPx, cenPy, 0} { Point{pFoilUpper[]}; }
pFoilLower[] = Symmetry {0, 1, 0, -cenPy} { Duplicata { Point{pFoilUpper[]}; } };
Delete{Point{pFoilLower[0]};}
Delete{Point{pFoilLower[nX]};}
pFoilLower[0] = pFoilUpper[0];
pFoilLower[nX] = pFoilUpper[nX];


lFoilUpper = newl;
Spline(lFoilUpper) = pFoilUpper[];
lFoilLower = newl;
Spline(lFoilLower) = pFoilLower[];

centerPoint = newp;
Point(centerPoint) = {cenPx,cenPy,0};

//pCircleEast  = newp; Point(pCircleEast)  = {cenPx+minor,cenPy,0};
pCircleWest  = newp; Point(pCircleWest)  = {x[0]-minor,y[0],0};
//pCircleNorth = newp; Point(pCircleNorth) = {cenPx,cenPy+major,0};
//pCircleSouth = newp; Point(pCircleSouth) = {cenPx,cenPy-major,0};
//
pNorthEast = newp; Point(pNorthEast) = {x[nX]+minor,y[nX]+major,0};
pSouthEast = newp; Point(pSouthEast) = {x[nX]+minor,y[nX]-major,0};
pEast = newp; Point(pEast) = {x[nX]+minor,y[nX],0};
//
pCircleNorth = newp; Point(pCircleNorth) = {x[nX],y[nX]+major,0};
pCircleSouth = newp; Point(pCircleSouth) = {x[nX],y[nX]-major,0};
//
//circleQ1 = newl; Ellipsis(circleQ1) = {pCircleEast , centerPoint, pCircleEast, pCircleNorth};
circleQ2 = newl; Ellipsis(circleQ2) = {pCircleNorth, pFoilUpper[nX], pCircleNorth, pCircleWest};
circleQ3 = newl; Ellipsis(circleQ3) = {pCircleWest , pFoilUpper[nX], pCircleWest, pCircleSouth};
//circleQ4 = newl; Ellipsis(circleQ4) = {pCircleSouth, centerPoint, pCircleSouth, pCircleEast};
//
lww = newl;
Line(lww) = {pCircleWest, pFoilUpper[0]};
//lnn = newl;
//Line(lnn) = {pCircleNorth, pFoilUpper[nX/2]};
//liwe = newl;
//For i In {0:nX}
	//temp[i] = pFoilUpper[i];
//EndFor
//Spline(liwe) = {temp[]};
//lee = newl;
//Line(lee) = {pCircleEast, pFoilUpper[nX]};
//line = newl;
//For i In {nX/2:nX}
	//temp[i-nX/2] = pFoilUpper[i];
//EndFor
//Spline(line) = {temp[]};
//lss = newl;
//Line(lss) = {pCircleSouth, pFoilLower[nX/2]};
//lise = newl;
//For i In {nX/2:nX}
//	temp[i-nX/2] = pFoilLower[i];
//EndFor
//Spline(lise) = {temp[]};
//liws = newl;
//For i In {0:nX/2}
//	temp[i] = pFoilLower[i];
//EndFor
//Spline(liws) = {temp[]};
//
len = newl;
Line(len) = {pFoilUpper[nX], pCircleNorth};
les = newl;
Line(les) = {pFoilUpper[nX], pCircleSouth};
lne = newl;
Line(lne) = {pCircleNorth, pNorthEast};
lfne = newl;
Line(lfne) = {pNorthEast, pEast};
lfee = newl;
Line(lfee) = {pEast, pFoilUpper[nX]};
lfse = newl;
Line(lfse) = {pSouthEast, pEast};
lse = newl;
Line(lse) = {pCircleSouth, pSouthEast};

llnw = newll;
Curve Loop(llnw) = {-lww, -circleQ2, -len, -lFoilUpper};
llne = newll;
Curve Loop(llne) = {len, lne, lfne, lfee};
llse = newll;
Curve Loop(llse) = {les, lse, lfse, lfee};
llsw = newll;
Curve Loop(llsw) = {les, -circleQ3, lww, lFoilLower};

Transfinite Line{len,les} = transY Using Progression 1/stretchVert;
Transfinite Line{lww,lfne,lfse,lfee} = transY Using Progression stretchVert;

Transfinite Line{lne,lse} = transWake Using Progression 1/stretchWake;
Transfinite Line{lfee} = transWake Using Progression stretchWake;

Transfinite Line{circleQ2, circleQ3, lFoilUpper,lFoilLower} = transX Using Bump 0.5;
//Transfinite Line{lFoilLower} = transX Using Bump 1/0.5;
//
Plane Surface(1) = {llnw};
Plane Surface(2) = {llne};
Plane Surface(3) = {llsw};
Plane Surface(4) = {llse};
//
Transfinite Surface{1};
Transfinite Surface{2};
Transfinite Surface{3};
Transfinite Surface{4};
//
Recombine Surface{1};
Recombine Surface{2};
Recombine Surface{3};
Recombine Surface{4};
//
out[] = Extrude {0,0,depth} { Surface{1,2,3,4}; Layers{nz}; Recombine;};
//
Printf("out[0], empty: %g",out[0]);
Printf("out[1], empty: %g",out[1]);
Printf("out[2], skip: %g",out[2]);
Printf("out[3], outer: %g",out[3]);
Printf("out[4], skip: %g",out[4]);
Printf("out[5], foil: %g",out[5]);

Printf("out[6], empty: %g",out[6]);
Printf("out[7], empty: %g",out[7]);
Printf("out[8], skip: %g",out[8]);
Printf("out[9], outer: %g",out[9]);
Printf("out[10], skip: %g",out[10]);
Printf("out[11], foil: %g",out[11]);

Printf("out[12], empty: %g",out[12]);
Printf("out[13], empty: %g",out[13]);
Printf("out[14], skip: %g",out[14]);
Printf("out[15], outer: %g",out[15]);
Printf("out[16], skip: %g",out[16]);
Printf("out[17], foil: %g",out[17]);

Printf("out[18], empty: %g",out[18]);
Printf("out[19], empty: %g",out[19]);
Printf("out[20], skip: %g",out[20]);
Printf("out[21], outer: %g",out[21]);
Printf("out[22], skip: %g",out[22]);
Printf("out[23], foil: %g",out[23]);
//
outerbc[0] = out[3];
outerbc[1] = out[9];
outerbc[2] = out[10];
outerbc[3] = out[15];
outerbc[4] = out[21];
outerbc[5] = out[22];

wallbc[0] = out[5];
wallbc[2] = out[17];

emptybc[0] = out[0];
emptybc[1] = out[1];
emptybc[2] = out[6];
emptybc[3] = out[7];
emptybc[4] = out[12];
emptybc[5] = out[13];
emptybc[6] = out[18];
emptybc[7] = out[19];

Printf("outerbc[0]: %g",outerbc[0]);
Printf("outerbc[1]: %g",outerbc[1]);
Printf("outerbc[2]: %g",outerbc[2]);
Printf("outerbc[3]: %g",outerbc[3]);
Printf("outerbc[4]: %g",outerbc[4]);
Printf("outerbc[5]: %g",outerbc[5]);

Printf("wallbc[0]: %g",wallbc[0]);
Printf("wallbc[1]: %g",wallbc[1]);

Printf("emptybc[0]: %g",emptybc[0]);
Printf("emptybc[1]: %g",emptybc[1]);
Printf("emptybc[2]: %g",emptybc[2]);
Printf("emptybc[3]: %g",emptybc[3]);
Printf("emptybc[4]: %g",emptybc[4]);
Printf("emptybc[5]: %g",emptybc[5]);

Physical Surface(1) = wallbc[];
Physical Surface(2) = outerbc[];
Physical Surface(3) = emptybc[];
Physical Volume(4) = {1,2,3,4};
