X = 0.5;
Y = 0.1;
Z = 0.1;
mx = 100;

p0 = newp; Point(p0) = {-X, -Y/2, 0};
p1 = newp; Point(p1) = {-X, Y/2, 0};
p2 = newp; Point(p2) = {-X, Y/2, Z};
p3 = newp; Point(p3) = {-X, -Y/2, Z};

l0 = newl; Line(l0) = {p0, p1};
l1 = newl; Line(l1) = {p1, p2};
l2 = newl; Line(l2) = {p2, p3};
l3 = newl; Line(l3) = {p3, p0};

Transfinite Line{l0, l1, l2, l3} = 1;

ll = newll; Line Loop(ll) = {l0, l1, l2, l3};
s0 = news; Plane Surface(s0) = {ll};

Transfinite Surface{s0};
Recombine Surface{s0};

out[] = Extrude {2*X,0,0} {Surface{s0}; Layers{mx}; Recombine;};

//Printf("out[0], : %g", out[0]);
//Printf("out[1], : %g", out[1]);
//Printf("out[2], : %g", out[2]);
//Printf("out[3], : %g", out[3]);
//Printf("out[4], : %g", out[4]);
//Printf("out[5], : %g", out[5]);

empty[] = {out[2], out[3], out[4], out[5]};
diric[] = {s0, out[0]};

Physical Surface(3) = {empty[]};
Physical Surface(2) = {diric[]};
Physical Volume(4) = {1};
