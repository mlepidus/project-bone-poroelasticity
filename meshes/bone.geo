SetFactory("OpenCASCADE");
Cylinder(1) = {1, 1, 0, 0, 0, 1, 0.7, 2*Pi};
Cylinder(2) = {1, 1, 0, 0, 0, 1, 1, 2*Pi};
s() = BooleanDifference{ Surface{2}; Delete; }{ Surface{1}; Delete; };
v() = BooleanDifference{ Volume{2}; Delete; }{ Volume{1}; Delete; };
MeshSize{ PointsOf{ Volume{2}; } } = 0.1;
Physical Volume("tutto",28) = {v(), 2};


//+
Physical Surface(29) = {3};
//+
Physical Surface(30) = {3};
//+
Physical Surface(31) = {3};
//+
Physical Surface(29) -= {3};
//+
Physical Surface(30) -= {3};
//+
Physical Surface(31) -= {3};
//+
Physical Surface(32) -= {3};
//+
Physical Surface(32) -= {3};
//+
Physical Surface(32) -= {2};
//+
Physical Surface(32) -= {4};
//+
Physical Surface(32) -= {1};
//+
Physical Surface("lower", 32) = {4};
//+
Physical Surface("outer", 33) = {2, 4};
//+
Physical Surface("inner", 34) = {1};
//+
Physical Surface("upper", 35) = {3};
//+
Physical Surface(" lower", 32) -= {3};
//+
Physical Surface(" lower", 32) -= {1};
//+
Physical Surface(" outer", 33) -= {3};
//+
Physical Surface(" lower", 32) -= {3};
//+
Physical Surface(" lower", 32) -= {3, 4};
//+
Physical Surface(" outer", 33) -= {4};
//+
Physical Surface("lower", 36) = {4};
