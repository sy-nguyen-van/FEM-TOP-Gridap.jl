// Parameters
length = 100;
height = 100;
cut_length = 60;
cut_height = 60;
length_F = 5;
thickness = 10; // Extrusion in z-direction

// Outer rectangle points (counter-clockwise in XY plane)
Point(1) = {0, 0, 0, 1.0};
Point(2) = {length, 0, 0, 1.0};
Point(3) = {length, height - cut_height - length_F, 0, 1.0};
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
Line Loop(1) = {1, 2, 3, 4, 5, 6, 7};
Plane Surface(1) = {1};

// Recombine for quad mesh
Mesh.Algorithm = 8;
Mesh.RecombineAll = 1;
Recombine Surface{1};

// Extrude to 3D solid
out[] = Extrude {0, 0, thickness} {
  Surface{1};
  Layers{thickness}; // Number of layers = thickness (adjust as needed)
  Recombine;
};

// Physical groups for boundary conditions
Physical Surface("Load") = {27}; // you can change to specific face of interest
Physical Surface("Top") = {39};  // you can change to specific face of interest

// Generate 3D mesh
Mesh 3;
Mesh.SaveAll = 0; // do not save lower-dim elements

Save "Lbracket3d_Sy.msh";