//SetFactory("OpenCASCADE");
Mesh.Algorithm = 8;  // Frontal Delaunay for quads
//Mesh.RecombinationAlgorithm = 2;

Point(1) = {0, 0, 0, 0.5};
Point(2) = {20, 0, 0, 0.5};
Point(3) = {20, 7, 0, 0.5};
Point(4) = {0, 7, 0, 0.5};
Point(5) = {1.5, 0, 0, 0.5};
Point(6) = {1.5, 1.5, 0, 0.5};
Point(7) = {0, 1.5, 0, 0.5};

Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 1};
Line(5) = {5, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 7};

Curve Loop(9) = {1, 2, 3, 4};
Curve Loop(10) = {1, 5, 6, 7, 8, 4};

Plane Surface(1) = {10, 9};
Plane Surface(2) = {9};

Recombine Surface{1};
Recombine Surface{2};

Physical Curve("HorEdges", 11) = {1, 5, 7};
Physical Curve("VerEdges", 12) = {2, 3, 4, 6, 8};
Physical Curve("Square", 13) = {1, 2, 3, 4};
Physical Surface("NoCorner", 3) = {1};
Physical Surface("Corner", 4) = {2}; 