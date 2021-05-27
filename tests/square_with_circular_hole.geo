lc = .5;

Point(1) = { 0, 0, 0, lc};
Point(2) = { 1, 0, 0, lc};
Point(3) = { 0, 1, 0, lc};
Point(4) = {-1, 0, 0, lc};
Point(5) = { 0,-1, 0, lc};
Point(6) = { 5,-5, 0, lc};
Point(7) = { 5, 5, 0, lc};
Point(8) = {-5, 5, 0, lc};
Point(9) = {-5,-5 ,0, lc};

// inner boundary
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};
Line Loop(1) = {1,2,3,4};
Physical Curve("Inner_Boundary") = {1,2,3,4};

// outer boundary
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {8,9};
Line(8) = {9,6};
Line Loop(2) = {5,6,7,8}; 
Physical Curve("Outer_Boundary_right") = {5};
Physical Curve("Outer_Boundary_top") = {6};
Physical Curve("Outer_Boundary_left") = {7};
Physical Curve("Outer_Boundary_bottom") = {8};

Plane Surface(1) = {2,1}; // line loop 2 - line loop 1
Physical Surface("SG") = {1};

