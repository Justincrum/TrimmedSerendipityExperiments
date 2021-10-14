
//Making a mesh for the monodomain cuboid.

SetFactory("OpenCASCADE");

//Mesh.Algorithm = 8;  // Frontal Delaunay for quads

// Inputs
fibreAxis = 20; //mm
NonFibreAxis = 7; //mm 
gridsize = 0.5;
meshThickness = 0.5; //mm

// Geometry of the corner cube
//Define the bottom
Point(1) = {0, 0, 0, gridsize};
Point(2) = {fibreAxis, 0, 0, gridsize};
Point(3) = {fibreAxis, NonFibreAxis, 0, gridsize};
Point(4) = {0, NonFibreAxis, 0, gridsize};
Point(5) = {1.5, 0, 0, gridsize};
Point(7) = {0, 1.5, 0, gridsize};
Point(6) = {1.5, 1.5, 0, gridsize};
Point(8) = {1.5, NonFibreAxis, 0, gridsize};

Line(1) = {1, 5};				// bottom left
Line(2) = {5, 2};				// bottom right
Line(3) = {2, 3};				// right
Line(4) = {3, 8};				// top right
Line(5) = {8, 4};               // top left
Line(6) = {4, 7};               // left top
Line(7) = {5, 6};               // corner right
Line(8) = {6, 7};               // corner top
Line(9) = {6, 8};               // mid top
Line(10) = {7, 1};              // left bottom
Line(11) = {8, 5};               // Middle line


Transfinite Curve {1} = 1.5/gridsize + 1;
Transfinite Curve {2} = (fibreAxis-1.5)/gridsize + 1;
Transfinite Curve {3} = 7/gridsize + 1;
Transfinite Curve {4} = (fibreAxis-1.5)/gridsize + 1;
Transfinite Curve {5} = 1.5/gridsize + 1;
Transfinite Curve {6} = (NonFibreAxis-1.5)/gridsize + 1;
Transfinite Curve {7} = 1.5/gridsize + 1;
Transfinite Curve {8} = 1.5/gridsize + 1;
Transfinite Curve {9} = 5.5/gridsize + 1;
Transfinite Curve {10} = 1.5/gridsize + 1;
Transfinite Curve {11} = 7/gridsize + 1;

Curve Loop(1) = {1, 7, 8, 10};
Plane Surface(1) = {1};
Transfinite Surface{1};
Recombine Surface{1};

Curve Loop(2) = {2, 3, 4, 11};
Plane Surface(2) = {2};
Transfinite Surface{2};
Recombine Surface{2};

Curve Loop(3) = {-8, 9, 5, 6};
Plane Surface(3) = {3};
Transfinite Surface{3};
Recombine Surface{3};
//Physical Surface(3) = {3};


surfaceVectorMesh[] = Extrude {0, 0, 3} {
Surface{2};
Layers{6};
Recombine;
};

//surfaceVectorMesh contains in the following order:
//[0]	- front surface (opposed to source surface)
//[1] - extruded volume
//[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
//[3] - right surface (belonging to 2nd line in "Line Loop (6)")
//[4] - top surface (belonging to 3rd line in "Line Loop (6)")
//[5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface("front", 4) = surfaceVectorMesh[0];
Physical Volume("internal", 1) = surfaceVectorMesh[1];
Physical Surface("bottom", 5) = surfaceVectorMesh[2];
Physical Surface("right", 6) = surfaceVectorMesh[3];
Physical Surface("top", 7) = surfaceVectorMesh[4];
Physical Surface("left", 8) = surfaceVectorMesh[5];
Physical Surface("back", 9) = {2};

surfaceVectorCorner[] = Extrude {0, 0, 1.5} {
Surface{1};
Layers{3};
Recombine;
};

//surfaceVectorMesh contains in the following order:
//[0]	- front surface (opposed to source surface)
//[1] - extruded volume
//[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
//[3] - right surface (belonging to 2nd line in "Line Loop (6)")
//[4] - top surface (belonging to 3rd line in "Line Loop (6)")
//[5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface("frontC", 10) = surfaceVectorCorner[0];
Physical Volume("internalC", 2) = surfaceVectorCorner[1];
Physical Surface("bottomC", 11) = surfaceVectorCorner[2];
Physical Surface("rightC", 12) = surfaceVectorCorner[3];
Physical Surface("topC", 13) = surfaceVectorCorner[4];
Physical Surface("leftC", 14) = surfaceVectorCorner[5];
Physical Surface("backC", 15) = {1};

surfaceVectorCornerTop[] = Extrude {0, 0, 1.5} {
Surface{13};
Layers{3};
Recombine;
};

//surfaceVectorMesh contains in the following order:
//[0]	- front surface (opposed to source surface)
//[1] - extruded volume
//[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
//[3] - right surface (belonging to 2nd line in "Line Loop (6)")
//[4] - top surface (belonging to 3rd line in "Line Loop (6)")
//[5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface("frontCT", 16) = surfaceVectorCornerTop[0];
Physical Volume("internalCT", 4) = surfaceVectorCornerTop[1];
Physical Surface("bottomCT", 17) = surfaceVectorCornerTop[2];
Physical Surface("rightCT", 18) = surfaceVectorCornerTop[3];
Physical Surface("topCT", 19) = surfaceVectorCornerTop[4];
Physical Surface("leftCT", 20) = surfaceVectorCornerTop[5];
Physical Surface("backCT", 21) = {13};

surfaceVectorTopleft[] = Extrude {0, 0, 3} {
Surface{3};
Layers{6};
Recombine;
};

//surfaceVectorMesh contains in the following order:
//[0]	- front surface (opposed to source surface)
//[1] - extruded volume
//[2] - bottom surface (belonging to 1st line in "Line Loop (6)")
//[3] - right surface (belonging to 2nd line in "Line Loop (6)")
//[4] - top surface (belonging to 3rd line in "Line Loop (6)")
//[5] - left surface (belonging to 4th line in "Line Loop (6)") 
Physical Surface("frontleft", 22) = surfaceVectorTopleft[0];
Physical Volume("internalleft", 3) = surfaceVectorTopleft[1];
Physical Surface("bottomleft", 23) = surfaceVectorTopleft[2];
Physical Surface("rightleft", 24) = surfaceVectorTopleft[3];
Physical Surface("topleft", 25) = surfaceVectorTopleft[4];
Physical Surface("leftleft", 26) = surfaceVectorTopleft[5];
Physical Surface("backleft", 27) = {3};