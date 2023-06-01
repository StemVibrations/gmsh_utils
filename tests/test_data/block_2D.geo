// Gmsh project: created with gmsh-3.0.6-Windows64
boxdim = 1;
column_height = 1;
gridsize = 1;

// Create 2D square mesh
Mesh.ElementOrder = 1;
Point(1) = {0, 0, 0, gridsize};
Point(2) = {boxdim, 0, 0, gridsize};
Point(3) = {boxdim, column_height, 0, gridsize};
Point(4) = {0, column_height, 0, gridsize};

// create lines
Line(5) = {1, 2};
Line(6) = {2, 3};
Line(7) = {3, 4};
Line(8) = {4, 1};

// create surface
Line Loop(9) = {5, 6, 7, 8};
Plane Surface(10) = 9;


Transfinite Surface{10};

// Define the physical groups
Physical Surface("group_1") = 10;
