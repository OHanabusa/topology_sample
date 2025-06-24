// Gmsh project created on Wed Apr 30 16:41:26 2025
SetFactory("OpenCASCADE");

L = 20.0 ;
H = 10.0 ;
H_press = 2.0 ;
Point(1) = {0.0, -0.5 * H, -0, 1.0};
Point(2) = {L, -0.5 * H, -0, 1.0};
Point(3) = {L, -0.5 * H_press, -0, 1.0};
Point(4) = {L, 0.5 * H_press, -0, 1.0};
Point(5) = {L, 0.5 * H, -0, 1.0};
Point(6) = {0.0, 0.5 * H, -0, 1.0};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

Curve Loop(1) = {6, 1, 2, 3, 4, 5};

Plane Surface(1) = {1};
Physical Curve(0) = {1, 2, 4, 5};
Physical Curve(1) = {6};
Physical Curve(2) = {3};
Physical Surface(0) = {1};
