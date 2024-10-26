#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <set>

using std::vector;

struct Mesh
{
   vector<double> rMain{ };
   vector<double> zMain{ };
   vector<vector<int>> areasMesh{ };
};

struct SplittingMesh
{
   vector<double> r{ };
   vector<double> z{ };
};

struct ShiftArrays
{
   vector<int> Ir{ };
   vector<int> Iz{ };
};

struct BoundaryConditions
{
   vector<vector<int>> conds{ };
};

struct SparseMatrix
{
   vector<int> ig{ }, jg{ };
   vector<double> ggl{ };
   vector<double> ggu{ };
   vector<double> di{ };
};

struct SLAE
{
   SparseMatrix A{ };
   vector<double> q{ };
   vector<double> b{ };
};

struct SLAEFull
{
   vector<vector<double>> A{ };
   vector<double> x{ };
   vector<double> b{ };
};

struct Vectors
{
   std::vector<double> x{};
   std::vector<double> b{};
   std::vector<double> r{};
   std::vector<double> s{};
   std::vector<double> r0{};
   std::vector<double> z{};
   std::vector<double> p{};
   std::vector<double> y{};
   std::vector<double> az{};
   std::vector<double> ay{};
};

struct MatrixM
{
   std::vector<double> di{};
   std::vector<double> au{};
   std::vector<double> al{};
};

struct Sources
{
   vector<vector<double>> A{ { 0, -500, 0 },
                             { 0, 0, 0 },
                             { 0, 500, 0 } };
   vector<vector<double>> B{ { 100, -500, 0 },
                             { 100, 0, 0 },
                             { 100, 500, 0 } };
};

struct Receivers
{
   vector<vector<double>> M{ { 200, 0, 0 },
                             { 500, 0, 0 },
                             { 1000, 0, 0 } };
   vector<vector<double>> N{ { 300, 0, 0 },
                             { 600, 0, 0 },
                             { 1100, 0, 0 } };
};

struct Matrices
{
   vector<vector<double>> G = { { 1, -1 },
                                { -1, 1 } };

   vector<vector<double>> M = { { 1. / 3., 1. / 6. },
                                { 1. / 6., 1. / 3. } };

   vector<vector<double>> H = { { 1. / 12., 1. / 12.},
                                { 1. / 12., 1. / 4. } };
};