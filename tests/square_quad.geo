Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Recombine Surface{1};

Physical Curve("unten") = {1};
Physical Curve("rechts") = {2};
Physical Curve("oben") = {3};
Physical Curve("links") = {4};
Physical Surface("Einheitsquadrat") = {1};

// command for second order serendipity elements:
// gmsh -2 square_quad.geo -setnumber Mesh.ElementOrder 2 -setnumber Mesh.SecondOrderIncomplete 1
