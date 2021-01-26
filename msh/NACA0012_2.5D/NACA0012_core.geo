c=1; // chord length
t=12 * c / 100; // thickness -- NACA0012
nX=1000; // number of 'x' points
depth=1; // extrusion length
transCircleEach=20;
transAfoilEach=50;
major=500;
minor=500;
originx=0;
originy=0;

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
	pointNum[i] = newp;
	Point( pointNum[i] ) = {x[i],y[i],0};
EndFor

For i In {0:nX}
	SplineList[i] = i+1;
EndFor

splineNum = newl;
Spline(splineNum) = SplineList[];

Translate {cenPx, cenPy, 0} { Line{splineNum}; }
symLine[] = Symmetry {0, 1, 0, -cenPy} { Duplicata { Line{splineNum}; } };

Transfinite Line{splineNum}  = transAfoilEach Using Bump 0.6;
Transfinite Line{symLine[0]} = transAfoilEach Using Bump 0.6;

afoilLoop = newll;
Line Loop (afoilLoop) = {splineNum, -symLine[0]};

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

Transfinite Line{circleQ1} = transCircleEach Using Progression 1.0;
Transfinite Line{circleQ2} = transCircleEach Using Progression 1.0;
Transfinite Line{circleQ3} = transCircleEach Using Progression 1.0;
Transfinite Line{circleQ4} = transCircleEach Using Progression 1.0;

circleLoop = newll;
Line Loop(circleLoop) = {circleQ1, circleQ2, circleQ3, circleQ4};
Plane Surface(0) = {circleLoop, afoilLoop};

out[]  = Extrude {0,0,depth} { Surface{0}; Layers{1}; Recombine;};

LatSurfCircle[0] = out[2];
LatSurfCircle[1] = out[3];
LatSurfCircle[2] = out[4];
LatSurfCircle[3] = out[5];

For i In {0:#out[]-7}
    SurfAirfoil[i] = out[i+6];
EndFor
