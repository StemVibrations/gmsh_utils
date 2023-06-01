// Gmsh project created on Thu Jun 01 14:45:20 2023
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.1};
//+
Point(2) = {2, 0, 0, 0.1};
//+
Point(3) = {2, 1, 0, 0.1};
//+
Point(4) = {1, 1, 0, 0.1};
//+
Point(5) = {0, 1, 0, 0.1};
//+
Point(6) = {2, 1.5, 0, 0.1};
//+
Point(7) = {1, 1.5, 0, 0.1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 1};
//+
Line(6) = {3, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 4};
//+
Curve Loop(1) = {1, 2, 3, 4, 5};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {3, -8, -7, -6};
//+
Surface(2) = {3};
