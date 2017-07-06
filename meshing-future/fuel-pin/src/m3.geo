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
Line(1) = {4, 17}; //Transfinite Line {1} = 10 Using Progression 1;
//+
Line(2) = {17, 3}; //Transfinite Line {2} = 10 Using Progression 1;
//+
Line(3) = {3, 14}; //Transfinite Line {3} = 10 Using Progression 1;
//+
Line(4) = {14, 2}; //Transfinite Line {4} = 10 Using Progression 1;
//+
Line(5) = {2, 16}; //Transfinite Line {5} = 10 Using Progression 1;
//+
Line(6) = {16, 1}; //Transfinite Line {6} = 10 Using Progression 1;
//+
Line(7) = {1, 15}; //Transfinite Line {7} = 10 Using Progression 1;
//+
Line(8) = {15, 4}; //Transfinite Line {8} = 10 Using Progression 1;
//+
Line(9) = {4, 9}; //Transfinite Line {9} = 10 Using Progression 1;
//+
Line(10) = {3, 8}; //Transfinite Line {10} = 10 Using Progression 1;
//+
Line(11) = {17, 11}; //Transfinite Line {11} = 10 Using Progression 1;
//+
Line(12) = {10, 16}; //Transfinite Line {12} = 10 Using Progression 1;
//+
Line(13) = {6, 1}; //Transfinite Line {13} = 10 Using Progression 1;
//+
Line(14) = {12, 15}; //Transfinite Line {14} = 10 Using Progression 1;
//+
Line(15) = {13, 14}; //Transfinite Line {15} = 10 Using Progression 1;
//+
Line(16) = {7, 2}; //Transfinite Line {16} = 10 Using Progression 1;
//+
Line(17) = {8, 13}; //Transfinite Line {17} = 10 Using Progression 1;
//+
Line(18) = {13, 7}; //Transfinite Line {18} = 10 Using Progression 1;
//+
Line(19) = {7, 10}; //Transfinite Line {19} = 10 Using Progression 1;
//+
Line(20) = {10, 6}; //Transfinite Line {20} = 10 Using Progression 1;
//+
Line(21) = {12, 6}; //Transfinite Line {21} = 10 Using Progression 1;
//+
Line(22) = {9, 12}; //Transfinite Line {22} = 10 Using Progression 1;
//+
Line(23) = {11, 9}; //Transfinite Line {23} = 10 Using Progression 1;
//+
Line(24) = {11, 8}; //Transfinite Line {24} = 10 Using Progression 1;
//+
Line(25) = {9, 5}; //Transfinite Line {25} = 10 Using Progression 1;
//+
Line(26) = {8, 5}; //Transfinite Line {26} = 10 Using Progression 1;
//+
Line(27) = {7, 5}; //Transfinite Line {27} = 10 Using Progression 1;
//+
Line(28) = {6, 5}; //Transfinite Line {28} = 10 Using Progression 1;
//+
Line Loop(29) = {24, -10, -2, 11};
//+
Plane Surface(30) = {29};
//+
Line Loop(31) = {10, 17, 15, -3};
//+
Plane Surface(32) = {31};
//+
Line Loop(33) = {15, 4, -16, -18};
//+
Plane Surface(34) = {33};
//+
Line Loop(35) = {5, -12, -19, 16};
//+
Plane Surface(36) = {35};
//+
Line Loop(37) = {6, -13, -20, 12};
//+
Plane Surface(38) = {37};
//+
Line Loop(39) = {13, 7, -14, 21};
//+
Plane Surface(40) = {39};
//+
Line Loop(41) = {8, 9, 22, 14};
//+
Plane Surface(42) = {41};
//+
Line Loop(43) = {1, 11, 23, -9};
//+
Plane Surface(44) = {43};
//+
Line Loop(45) = {23, 25, -26, -24};
//+
Plane Surface(46) = {45};
//+
Line Loop(47) = {26, -27, -18, -17};
//+
Plane Surface(48) = {47};
//+
Line Loop(49) = {27, -28, -20, -19};
//+
Plane Surface(50) = {49};
//+
Line Loop(51) = {28, -25, 22, 21};
//+
Plane Surface(52) = {51};
//+
//Transfinite Surface {30};
//+
//Transfinite Surface {32};
//+
//Transfinite Surface {34};
//+
//Transfinite Surface {36};
//+
//Transfinite Surface {38};
//+
//Transfinite Surface {40};
//+
//Transfinite Surface {42};
//+
//Transfinite Surface {44};
//+
//Transfinite Surface {46};
//+
//Transfinite Surface {48};
//+
//Transfinite Surface {50};
//+
//Transfinite Surface {52};
//+
Recombine Surface {30};
//+
Recombine Surface {32};
//+
Recombine Surface {34};
//+
Recombine Surface {36};
//+
Recombine Surface {38};
//+
Recombine Surface {40};
//+
Recombine Surface {42};
//+
Recombine Surface {44};
//+
Recombine Surface {46};
//+
Recombine Surface {48};
//+
Recombine Surface {50};
//+
Recombine Surface {52};

Mesh.CharacteristicLengthMin = 7.2;
//Mesh.CharacteristicLengthFactor = 9;
Mesh.SubdivisionAlgorithm = 1;
Mesh.Smoothing = 20;
//Mesh.RecombineAll = 1;
