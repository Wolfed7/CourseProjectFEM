using CourseProjectFEM;
using System.Globalization;

CultureInfo.CurrentCulture = new CultureInfo("en-US");

const string path1 = "Mesh/Input/Settings.dat";
const string path2 = "Mesh/Input/Splitting.dat";
const string path3 = "Mesh/Input/BoundaryConditions.dat";
const string path5 = "Mesh/Input/Time.dat";

const string path4 = "Mesh/Points.dat";

var mesh = new Mesh();
mesh.Input(path1, path2, path3, path5);
mesh.Build();
mesh.Output(path4);

var fem = new FEM();
fem.SetMesh(mesh);
fem.SetSolver(new BCG(1e-15, 1500));

fem.Compute();