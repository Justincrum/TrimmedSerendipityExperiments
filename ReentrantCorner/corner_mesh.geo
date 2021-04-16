//Create all the vertices.
Point(1) = {-1, 1, 0, 0.05};
Point(2) = {-1, 0, 0, 0.5};
Point(3) = {-1, -1, 0, 0.5};
Point(4) = {0, -1, 0, 0.5};
Point(5) = {0,  0, 0, 0.05};
Point(6) = {1,  0, 0, 0.25};
Point(7) = {1, 1, 0, 0.5};
Point(8) = {0, 1, 0, 0.05};

//Create all the edges.
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 1};

//Define the curves.
Curve Loop(7) = {1, 2, 3, 4, 5, 6, 7, 8};

Plane Surface(1) = {7};

Physical Curve("ZeroBCs", 8) = {4, 5};
Physical Curve("YNegative", 9) = {2, 3};
Physical Curve("YPositive", 10) = {1, 6, 7, 8};
Physical Surface("ReentrantCorner", 13) = {1};

Recombine Surface{:};

Mesh.Algorithm = 8;  // Frontal Delaunay for quads
Mesh.RecombinationAlgorithm = 2;
// RecombineMesh;

