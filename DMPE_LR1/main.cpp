#include "functions.h"
#include "iomanip"

int main()
{
   Mesh mesh{ };
   SplittingMesh sMesh{ };
   ShiftArrays shifts{ };
   BoundaryConditions conditions{ };
   SLAE slae{ };
   FunctionProblem funcs{ };
   MatrixM M{ };
   Vectors vectors{ };

   vector<double> I0 = { 4.5, 5.5, 7.5 };
   vector<double> Iresult;

   inputMesh(mesh);

   std::ifstream streamSplitMesh("splittingMesh.txt");
   inputSplittingMesh(mesh.rMain, sMesh.r, shifts.Ir, streamSplitMesh);
   inputSplittingMesh(mesh.zMain, sMesh.z, shifts.Iz, streamSplitMesh);
   streamSplitMesh.close();

   inputBoundaryCondtitions(conditions);

   generatePortraitOfSparseMatrix(slae, sMesh.r.size(), sMesh.z.size());

   calcGlobalMatrix(slae, sMesh, shifts, mesh, funcs, conditions);
//   LOS_LU(slae, 10000, 1e-15);
//   GMRES(slae, 10, 1000, 1e-14);
   BiCGStabWithLUPrecond(slae, M, vectors, 1e-15, 10000);
//   std::cout << std::setprecision(15) << NumericalValueAtPoint(sMesh, slae.q, 500, 0);

   solutionInverseProblem(sMesh, slae.q, funcs, I0, Iresult);

   return 0;
}