Point(1) = {1, 0, 0, 1.0};
Point(2) = {1, 1, 0, 1.0};
Point(3) = {0, 1, 0, 1.0};
Point(4) = {0, 0, 1, 1.0};
Point(5) = {1, 0, 1, 1.0};
Point(6) = {1, 1, 1, 1.0};
Point(7) = {0, 1, 1, 1.0};
Point(8) = {0, 0, 0, 1.0};

Line(1) = {7, 6};
Line(2) = {6, 5};
Line(3) = {5, 1};
Line(4) = {1, 8};
Line(5) = {8, 3};
Line(6) = {3, 7};
Line(7) = {7, 4};
Line(8) = {4, 8};
Line(9) = {4, 5};
Line(10) = {2, 1};
Line(11) = {2, 6};
Line(12) = {2, 3};

Line Loop(1) = {6, 1, -11, 12};
Plane Surface(1) = {1};
Line Loop(2) = {11, 2, 3, -10};
Plane Surface(2) = {2};
Line Loop(3) = {2, -9, -7, 1};
Plane Surface(3) = {-3};
Line Loop(4) = {6, 7, 8, 5};
Plane Surface(4) = {-4};
Line Loop(5) = {8, -4, -3, -9};
Plane Surface(5) = {5};
Line Loop(6) = {10, 4, 5, -12};
Plane Surface(6) = {6};
Surface Loop(1) = {6, 2, 1, 4, 3, 5};
Volume(1) = {1};

Physical Surface("vier") = {4, 3, 2, 6};
Physical Surface("A") = {1};
Physical Surface("B") = {5};
Physical Volume("Wuerfel") = {1};

// Lay1 = gmsh.model.addPhysicalGroup(dim, [1])
// gmsh.model.setPhysicalName(dim, Lay1, "SedimentLayer1")
