// Gmsh project created on Thu Jun 22 04:22:11 2017
spacing = 5.0;
//+
Point(1) = {1, 1, 0, 1};
//+
Point(2) = {1, -1, 0, 1};
//+
Point(3) = {-1, -1, 0, 1};
//+
Point(4) = {-1, 1, 0, 1};
//+
Point(5) = {0, 0, 0, 1};
//+
Point(6) = {0.56568542494, 0.56568542494, 0, spacing};
//+
Point(7) = {0.56568542494, -0.56568542494, 0, spacing};
//+
Point(8) = {-0.56568542494, -0.56568542494, 0, spacing};
//+
Point(9) = {-0.56568542494, 0.56568542494, 0, spacing};
//+
Point(10) = {0.8, 0, 0, spacing};
//+
Point(11) = {-0.8, 0, 0, spacing};
//+
Point(12) = {-0., 0.8, 0, spacing};
//+
Point(13) = {-0., -0.8, 0, spacing};
//+
Point(14) = {-0., -1, 0, spacing};
//+
Point(15) = {-0., 1, 0, spacing};
//+
Point(16) = {1, 0, 0, spacing};
//+
Point(17) = {-1, 0, 0, spacing};
//+
Line(1) = {4, 17};
//+
Line(2) = {17, 3};
//+
Line(3) = {3, 14};
//+
Line(4) = {14, 2};
//+
Line(5) = {2, 16};
//+
Line(6) = {16, 1};
//+
Line(7) = {1, 15};
//+
Line(8) = {15, 4};
//+
Circle(9) = {9, 5, 11};
//+
Circle(10) = {11, 5, 8};
//+
Circle(11) = {8, 5, 13};
//+
Circle(12) = {13, 5, 7};
//+
Circle(13) = {7, 5, 10};
//+
Circle(14) = {10, 5, 6};
//+
Circle(15) = {6, 5, 12};
//+
Circle(16) = {12, 5, 9};
//+
Line(17) = {4, 9};
//+
Line(18) = {1, 6};
//+
Line(19) = {7, 2};
//+
Line(20) = {8, 3};
//+
Line Loop(21) = {17, 9, 10, 20, -2, -1};
//+
Plane Surface(22) = {21};
//+
Line Loop(23) = {20, 3, 4, -19, -12, -11};
//+
Plane Surface(24) = {23};
//+
Line Loop(25) = {19, 5, 6, 18, -14, -13};
//+
Plane Surface(26) = {25};
//+
Line Loop(27) = {18, 15, 16, -17, -8, -7};
//+
Plane Surface(28) = {27};
//+
Line Loop(29) = {9, 10, 11, 12, 13, 14, 15, 16};
//+
Plane Surface(30) = {29};
//+
Recombine Surface {22};
//+
Recombine Surface {24};
//+
Recombine Surface {26};
//+
Recombine Surface {28};
//+
Recombine Surface {30};
Mesh.ElementOrder = 1;
// Mesh.RecombineAll = 1;
