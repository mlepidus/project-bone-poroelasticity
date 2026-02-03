Mesh.MshFileVersion = 2.2;
SetFactory("OpenCASCADE");

// -----------------------------
// Parametri geometrici
// -----------------------------
xc = 1;
yc = 1;
z0 = 0;
h  = 1;

r_in  = 0.7;
r_out = 1.0;

// -----------------------------
// Solidi pieni
// -----------------------------
Cylinder(1) = {xc, yc, z0, 0, 0, h, r_out, 2*Pi};
Cylinder(2) = {xc, yc, z0, 0, 0, h, r_in,  2*Pi};

// -----------------------------
// Cilindro cavo (volume unico)
// -----------------------------
hollow[] = BooleanDifference{ Volume{1}; Delete; }
                           { Volume{2}; Delete; };

// -----------------------------
// Mesh size
// -----------------------------
MeshSize{ PointsOf{ Volume{hollow[]}; } } = 0.1;

// -----------------------------
// Physical Volume
// -----------------------------
Physical Volume("solid") = {hollow[]};

// -----------------------------
// Physical Surfaces
// -----------------------------
// Recupero automatico delle superfici del volume
surfs[] = Boundary{ Volume{hollow[]}; };

// Identificazione:
// - z = z0        -> lower
// - z = z0 + h    -> upper
// - r = r_out     -> outer
// - r = r_in      -> inner

Physical Surface("lower") = { Surface In BoundingBox
  {xc-r_out-1e-6, yc-r_out-1e-6, z0-1e-6,
   xc+r_out+1e-6, yc+r_out+1e-6, z0+1e-6} };

Physical Surface("upper") = { Surface In BoundingBox
  {xc-r_out-1e-6, yc-r_out-1e-6, z0+h-1e-6,
   xc+r_out+1e-6, yc+r_out+1e-6, z0+h+1e-6} };

Physical Surface("outer") = { Surface In BoundingBox
  {xc-r_out-1e-6, yc-r_out-1e-6, z0-1e-6,
   xc+r_out+1e-6, yc+r_out+1e-6, z0+h+1e-6} };

Physical Surface("inner") = { Surface In BoundingBox
  {xc-r_in-1e-6, yc-r_in-1e-6, z0-1e-6,
   xc+r_in+1e-6, yc+r_in+1e-6, z0+h+1e-6} };
