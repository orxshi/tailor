SetFactory("OpenCASCADE");

c=1; // chord length
t=12 * c / 100; // thickness -- NACA0012
nX=1000; // number of 'x' points
depth=10; // extrusion length
transCircleEach=12;
transVert=20;
major=10;
minor=10;
originx=0;
originy=0;
zoff=5;
nz = 1;

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

pCircleEast  = newp; Point(pCircleEast)  = {cenPx+minor,cenPy,0};
pCircleWest  = newp; Point(pCircleWest)  = {cenPx-minor,cenPy,0};
pCircleNorth = newp; Point(pCircleNorth) = {cenPx,cenPy+major,0};
pCircleSouth = newp; Point(pCircleSouth) = {cenPx,cenPy-major,0};

circleQ1 = newl; Ellipsis(circleQ1) = {pCircleEast , centerPoint, pCircleEast, pCircleNorth};
circleQ2 = newl; Ellipsis(circleQ2) = {pCircleNorth, centerPoint, pCircleNorth, pCircleWest};
circleQ3 = newl; Ellipsis(circleQ3) = {pCircleWest , centerPoint, pCircleWest, pCircleSouth};
circleQ4 = newl; Ellipsis(circleQ4) = {pCircleSouth, centerPoint, pCircleSouth, pCircleEast};

lww = newl;
Line(lww) = {pCircleWest, pFoilUpper[0]};
lnn = newl;
Line(lnn) = {pCircleNorth, pFoilUpper[nX/2]};
liwn = newl;
For i In {0:nX/2}
	temp[i] = pFoilUpper[i];
EndFor
Spline(liwn) = {temp[]};
lee = newl;
Line(lee) = {pCircleEast, pFoilUpper[nX]};
line = newl;
For i In {nX/2:nX}
	temp[i-nX/2] = pFoilUpper[i];
EndFor
Spline(line) = {temp[]};
lss = newl;
Line(lss) = {pCircleSouth, pFoilLower[nX/2]};
lise = newl;
For i In {nX/2:nX}
	temp[i-nX/2] = pFoilLower[i];
EndFor
Spline(lise) = {temp[]};
liws = newl;
For i In {0:nX/2}
	temp[i] = pFoilLower[i];
EndFor
Spline(liws) = {temp[]};

llnw = newll;
Curve Loop(llnw) = {-lww, -circleQ2, lnn, -liwn};
llne = newll;
Curve Loop(llne) = {-lnn, -circleQ1, lee, -line};
llse = newll;
Curve Loop(llse) = {-lee, -circleQ4, lss, lise};
llsw = newll;
Curve Loop(llsw) = {-lss, -circleQ3, lww, liws};

Plane Surface(1) = {llnw};
Plane Surface(2) = {llne};
Plane Surface(3) = {llsw};
Plane Surface(4) = {llse};

out[] = Extrude {0,0,depth} { Surface{1,2,3,4};};

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


loop_foil = newll;
//Curve Loop(loop_foil) = {21, 43, 87, 65};
Curve Loop(loop_foil) = {9, 11, 14, 13};

sur_foil = news;
Plane Surface(sur_foil) = {loop_foil};

out4[] = Extrude {0,0,zoff} {Surface{sur_foil};};

Printf("out4[0], sur_foil: %g",sur_foil);
Printf("out4[0], out4[0]: %g",out4[0]);
Printf("out4[1], out4[1]: %g",out4[1]);
Printf("out4[2], out4[2]: %g",out4[2]);
Printf("out4[3], out4[3]: %g",out4[3]);
Printf("out4[4], out4[4]: %g",out4[4]);
Printf("out4[5], out4[5]: %g",out4[5]);

Mesh.CharacteristicLengthExtendFromBoundary = 0;
Mesh.CharacteristicLengthFromPoints = 0;
Mesh.CharacteristicLengthFromCurvature = 0;

lc = 100;

foil_center = newp;
Point(foil_center) = {cenPx+c/2,cenPy,0};

Field[1] = Attractor;
//Field[1].NodesList = {foil_center};
Field[1].FacesList = {39, 61, 83, 105};
Field[1].NNodesByEdge = {100};

Field[2] = MathEval;
Field[2].F = Sprintf("F1/4 + %g", lc / 1000);

Background Field = 2;

cl_wh = newll;
Curve Loop(cl_wh) = {circleQ1, circleQ2, circleQ3, circleQ4};

sur_wh = news;
Plane Surface(sur_wh) = {cl_wh};

out6[] = Extrude {0,0,-zoff} {Surface{sur_wh};};

Printf("out6[0], sur_wh: %g",sur_wh);
Printf("out6[0], out4[0]: %g",out6[0]);
Printf("out6[1], out4[1]: %g",out6[1]);
Printf("out6[2], out4[2]: %g",out6[2]);
Printf("out6[3], out4[3]: %g",out6[3]);
Printf("out6[4], out4[4]: %g",out6[4]);
Printf("out6[5], out4[5]: %g",out6[5]);

wallbc[0] = out[5];
wallbc[1] = out[11];
wallbc[2] = out[17];
wallbc[3] = out[23];

emptybc[0] = out4[0];
emptybc[1] = out6[0];

outerbc[0] = out4[2];
outerbc[1] = out4[3];
outerbc[2] = out4[4];
outerbc[3] = out4[5];
outerbc[4] = out6[2];
outerbc[5] = out6[3];
outerbc[6] = out6[4];
outerbc[7] = out6[5];
outerbc[9] = out[3];
outerbc[10] = out[9];
outerbc[11] = out[15];
outerbc[12] = out[21];

Compound Surface{out[0], out[6], out[12], out[18], sur_foil};
Compound Surface{out[1], out[7], out[13], out[19], sur_wh};
