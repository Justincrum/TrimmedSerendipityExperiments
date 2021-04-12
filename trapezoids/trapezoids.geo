//Create all the vertices.
Point(1) = {0, 0, 0, 0.37};
Point(2) = {1, 0, 0, 0.37};
Point(3) = { 1, 1, 0, 0.37};
Point(4) = { 0,  1, 0, 0.37};
Point(5) = { 0,  0.15, 0, 0.37};
Point(6) = {0, 0.5, 0, 0.37};
Point(7) = {0.25, 0.5, 0, 0.37};
Point(8) = {0.25, 0.35, 0, 0.37};
Point(9) = {0.25, 0, 0, 0.37};
Point(10) = {0.5, 0, 0, 0.37};
Point(11) = {0.75, 0, 0, 0.37};
Point(12) = {0.5, 0.5, 0, 0.37};
Point(13) = {0.75, 0.5, 0, 0.37};
Point(14) = {1, 0.5, 0, 0.37};
Point(15) = {0.5, 0.15, 0, 0.37};
Point(16) = {0.75, 0.35, 0, 0.37};
Point(17) = {1, 0.15, 0, 0.37};
Point(18) = {0, 0.85, 0, 0.37};
Point(19) = {0.25, 0.65, 0, 0.37};
Point(20) = {0.5, 0.85, 0, 0.37};
Point(21) = {0.75, 0.65, 0, 0.37};
Point(22) = {1, 0.85, 0, 0.37};
Point(23) = {0.25, 1, 0, 0.37};
Point(24) = {0.5, 1, 0, 0.37};
Point(25) = {0.75, 1, 0, 0.37};

//Create all the edges.
Line(1) = {1, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 8};
Line(5) = {8, 5};
Line(6) = {8, 9};
Line(7) = {9, 1};
Line(8) = {9, 10};
Line(9) = {10, 11};
Line(10) = {11, 2};
Line(11) = {7, 12};
Line(12) = {12, 13};
Line(13) = {13, 14};
Line(14) = {8, 15};
Line(15) = {15, 16};
Line(16) = {16, 17};
Line(17) = {15, 10};
Line(18) = {12, 15};
Line(19) = {13, 16};
Line(20) = {16, 11};
Line(21) = {14, 17};
Line(22) = {17, 2};
Line(23) = {6, 18};
Line(24) = {18, 4};
Line(25) = {18, 19};
Line(26) = {19, 20};
Line(27) = {20, 21};
Line(28) = {21, 22};
Line(29) = {7, 19};
Line(30) = {12, 20};
Line(31) = {13, 21};
Line(32) = {14, 22};
Line(33) = {4, 23};
Line(34) = {23, 24};
Line(35) = {24, 25};
Line(36) = {25, 3};
Line(37) = {3, 22};
Line(38) = {23, 19};
Line(39) = {24, 20};
Line(40) = {21, 25};

//Define the curves.
Curve Loop(1) = {-5, 6, 7, 1};
Curve Loop(2) = {3, 4, 5, 2};
Curve Loop(3) = {25, -29, -3, 23};
Curve Loop(4) = {33, 38, -25, 24};
Curve Loop(5) = {34, 39, -26, -38};
Curve Loop(6) = {26, -30, -11, 29};
Curve Loop(7) = {11, 18, -14, -4};
Curve Loop(8) = {14, 17, -8, -6};
Curve Loop(9) = {15, 20, -9, -17};
Curve Loop(10) = {12, 19, -15, -18};
Curve Loop(11) = {27, -31, -12, 30};
Curve Loop(12) = {35, -40, -27, -39};
Curve Loop(13) = {36, 37, -28, 40};
Curve Loop(14) = {28, -32, -13, 31};
Curve Loop(15) = {13, 21, -16, -19};
Curve Loop(16) = {16, 22, -10, -20};

//Define the trapezoids based off the curve loops.
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};
Plane Surface(8) = {8};
Plane Surface(9) = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};
Plane Surface(12) = {12};
Plane Surface(13) = {13};
Plane Surface(14) = {14};
Plane Surface(15) = {15};
Plane Surface(16) = {16};

Recombine Surface{:};

//Physical Curve("HorEdges", 11) = {1, 3};
//Physical Curve("VerEdges", 12) = {2, 4};
//Physical Curve("Circle", 13) = {8, 7, 6, 5};
//Physical Surface("PunchedDom", 3) = {1};
//Physical Surface("Disc", 4) = {2};

Physical Surface("Trap1", 1) = {1};
Physical Surface("Trap2", 2) = {2};
Physical Surface("Trap3", 3) = {3};
Physical Surface("Trap4", 4) = {4};
Physical Surface("Trap5", 5) = {5};
Physical Surface("Trap6", 6) = {6};
Physical Surface("Trap7", 7) = {7};
Physical Surface("Trap8", 8) = {8};
Physical Surface("Trap9", 9) = {9};
Physical Surface("Trap10", 10) = {10};
Physical Surface("Trap11", 11) = {11};
Physical Surface("Trap12", 12) = {12};
Physical Surface("Trap13", 13) = {13};
Physical Surface("Trap14", 14) = {14};
Physical Surface("Trap15", 15) = {15};
Physical Surface("Trap16", 16) = {16};

//Recombine Surface{:};

Mesh.Algorithm = 8;  // Frontal Delaunay for quads
Mesh.RecombinationAlgorithm = 2;
// RecombineMesh;

