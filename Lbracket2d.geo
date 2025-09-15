// Parameters
length = 100;
height = 100;
cut_length = 60;
cut_height = 60;
length_F = 5;
// Outer rectangle points (counter-clockwise)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {length, 0, 0, 1.0};
Point(3) = {length, height - cut_height- length_F, 0, 1.0};
Point(4) = {length, height - cut_height , 0, 1.0};
Point(5) = {length - cut_length, height - cut_height, 0, 1.0};
Point(6) = {length - cut_length, height, 0, 1.0};
Point(7) = {0, height, 0, 1.0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 1};
// Line loop and surface
Line Loop(1) = {1, 2, 3, 4, 5, 6,7};
Plane Surface(1) = {1};

Physical Curve("Load") = {3};
Physical Curve("Top") = {6};

// Mesh settings
Mesh.Algorithm = 8; // Frontal-Delaunay for Quads
Mesh.RecombineAll = 1; // Enable quad meshing
Recombine Surface{1}; // Explicitly recombine the surface

// Generate mesh
Mesh 2;
