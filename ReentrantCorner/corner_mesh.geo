//Create all the vertices.
Point(1) = {-1, 1, 0, 0.05};
Point(2) = {-1, -1, 0, 0.5};
Point(3) = {0, -1, 0, 0.5};
Point(4) = {0,  0, 0, 0.05};
Point(5) = {1,  0, 0, 0.5};
Point(6) = {1, 1, 0, 0.5};

//Create all the edges.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

//Define the curves.
Curve Loop(7) = {1, 2, 3, 4, 5, 6};

Plane Surface(1) = {7};

Physical Curve("ZeroBCs", 8) = {3, 4};
Physical Curve("FunctionBCs", 9) = {1, 2, 5, 6};
Physical Surface("ReentrantCorner", 13) = {1};

Recombine Surface{:};

Mesh.Algorithm = 8;  // Frontal Delaunay for quads
Mesh.RecombinationAlgorithm = 2;
// RecombineMesh;

