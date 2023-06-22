// Gmsh project created on Thu Jun 22 14:24:42 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {4, 0, 0, 0.1};
//+
Point(3) = {4, 2, 0, 0.1};
//+
Point(4) = {0, 2, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Curve Loop(1) = {3, 4, 1, 2};
//+
Plane Surface(1) = {1};
//+
Point(5) = {1, 2, 0, 0.1};
//+
Point(6) = {3, 2, 0, 0.1};
//+
Point(7) = {3, 3, 0, 0.1};
//+
Point(8) = {1, 3, 0, 0.1};
//+
Line(5) = {5, 6};
//+
Line(6) = {6, 7};
//+
Line(7) = {7, 8};
//+
Line(8) = {8, 5};
//+
Curve Loop(2) = {7, 8, 5, 6};
//+
Plane Surface(2) = {2};
