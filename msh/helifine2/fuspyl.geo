Merge "fuspyl.brep";

Physical Surface(9) = {1};
Physical Surface(1) = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 4, 5, 6, 7, 8, 9};
Physical Volume(4) = {1};

Mesh.ElementOrder = 1;
Mesh.Algorithm = 1;
Mesh.Algorithm3D = 10;
Mesh.Format = 1;
Mesh.MshFileVersion = 2.2;
Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.RandomFactor = 1e-6;

lc = 0.006;
Field[1] = Distance;
Field[1].SurfacesList = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 4, 5, 6, 7, 8, 9};
Field[1].NumPointsPerCurve = 200;
Field[2] = MathEval;
Field[2].F = Sprintf("F1/5 + %g", lc);
Background Field = 2;
