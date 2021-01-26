// geo file for meshing with Gmsh meshing software created by FreeCAD

// open brep geometry
Merge "hubshaft.brep";

// group data
Physical Surface(11) = {1};
Physical Volume(4) = {1};
Physical Surface(1) = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 4, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 6, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 7, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 8, 80, 9};

// Characteristic Length
// no boundary layer settings for this mesh
// min, max Characteristic Length
Mesh.CharacteristicLengthMax = 1e+22;
Mesh.CharacteristicLengthMin = 0.0;

// optimize the mesh
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 0;
Mesh.ElementOrder = 1;

// mesh order
Mesh.ElementOrder = 1;

// mesh algorithm, only a few algorithms are usable with 3D boundary layer generation
// 2D mesh algorithm (1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=BAMG, 8=DelQuad)
Mesh.Algorithm = 2;
// 3D mesh algorithm (1=Delaunay, 2=New Delaunay, 4=Frontal, 5=Frontal Delaunay, 6=Frontal Hex, 7=MMG3D, 9=R-tree)
Mesh.Algorithm3D = 1;

// meshing
Geometry.Tolerance = 1e-06; // set geometrical tolerance (also used for merging nodes)

// save
Mesh.Format = 1;
// Ignore Physical definitions and save all elements;


//////////////////////////////////////////////////////////////////////
// Gmsh documentation:
// http://gmsh.info/doc/texinfo/gmsh.html#Mesh
//
// We do not check if something went wrong, like negative jacobians etc. You can run Gmsh manually yourself: 
//
// to see full Gmsh log, run in bash:
// /home/orhan/gmsh-4.0.1-Linux64/bin/gmsh - /tmp/shape2mesh.geo
//
// to run Gmsh and keep file in Gmsh GUI (with log), run in bash:
// /home/orhan/gmsh-4.0.1-Linux64/bin/gmsh /tmp/shape2mesh.geo
Mesh.MshFileVersion = 2.2;

lc = 1;
Field[1] = Distance;
Field[1].FacesList = {10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 2, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 3, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 4, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 5, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 6, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 7, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 8, 80, 9};
Field[1].NNodesByEdge = 100;
Field[2] = MathEval;
Field[2].F = Sprintf("F1/2 + %g", lc / 1000);
Background Field = 2;
