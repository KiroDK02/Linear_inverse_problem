#include "functions.h"

void calcGlobalMatrix(SLAE &slae, SplittingMesh &sMesh, ShiftArrays &shifts, Mesh &mesh, FunctionProblem &funcs, BoundaryConditions &conditions)
{
   const double PI = 3.14159265358979;

   auto &f = funcs.f;
   auto &r = sMesh.r;
   auto &z = sMesh.z;

   const int rSize = r.size();
   const int zSize = z.size();

   int l = 0;

   for (int s = 0; s < zSize - 1; s++)
   {
      for (int p = 0; p < rSize - 1; p++)
      {
         int numArea = numberOfCurrentArea(shifts, mesh.areasMesh, p, s, l);
         
         vector<int> globalNumbers{ rSize * s + p, rSize * s + p + 1,
                                    rSize * (s + 1) + p, rSize * (s + 1) + p + 1 };

         double hr = r[p + 1] - r[p];
         double hz = z[s + 1] - z[s];

         if (numArea != -1)
         {
            vector<vector<double>> stiffness{ };
            vector<double> localRightPart(4);

            setSizeMatrix(4, 4, stiffness);

 //           vector<double> valueF{ f[numArea](r[p], z[s]), f[numArea](r[p + 1], z[s]),
 //                                  f[numArea](r[p], z[s + 1]), f[numArea](r[p + 1], z[s + 1]) };

            calcLocalStiffnessMatrix(r[p], hr, hz, stiffness);
 //           calcLocalRightPart(valueF, r[p], hr, hz, localRightPart);

            addLocalMatrixToGlobalMatrix(globalNumbers, funcs.sigma[numArea], stiffness, slae);

 //           addLocalRightPartToGlobal(globalNumbers, localRightPart, slae);
         }
         else
         {
            vector<vector<int>> localNums = { { p, s }, { p + 1, s }, { p, s + 1 }, { p + 1, s + 1 } };

            int l1 = 0;
            for (int i = 0; i < 4; i++)
               if (isFictitiousNode(shifts, mesh.areasMesh, localNums[i][0], localNums[i][1], l1))
                  slae.b[globalNumbers[i]] = 1;
         }
      }
   }

   slae.b[rSize * (zSize - 1)] = 1 / (2 * PI);

   addSecondBoundaryConditions(shifts, sMesh, slae, conditions, funcs);
   addFirstBoundaryConditions(shifts, sMesh, slae, conditions, funcs);
}

void addFirstBoundaryConditions(ShiftArrays &shifts, SplittingMesh &sMesh, SLAE &slae, BoundaryConditions &conditions, FunctionProblem &funcs)
{
   auto &conds = conditions.conds;
   auto &ug = funcs.ug;

   const int rSize = sMesh.r.size();
   const int zSize = sMesh.z.size();
   const int countConds = conds.size();

   for (int k = 0; k < countConds; k++)
   {
      int typeCondition = conds[k][0];

      if (typeCondition == 1)
      {
         int numFunction = conds[k][1];
         int p0 = conds[k][2];
         int p1 = conds[k][3];
         int s0 = conds[k][4];
         int s1 = conds[k][5];

         if (p0 == p1)
         {
            int p = shifts.Ir[p0];

            s0 = shifts.Iz[s0];
            s1 = shifts.Iz[s1];

            for (int s = s0; s <= s1; s++)
            {
               int globalNum = rSize * s + p;

               setStringMatrixInZero(slae.A, globalNum);
               slae.b[globalNum] = ug[numFunction](sMesh.r[p], sMesh.z[s]);
            }
         }
         else
         {
            int s = shifts.Iz[s0];

            p0 = shifts.Ir[p0];
            p1 = shifts.Ir[p1];

            for (int p = p0; p <= p1; p++)
            {
               int globalNum = rSize * s + p;

               setStringMatrixInZero(slae.A, globalNum);
               slae.b[globalNum] = ug[numFunction](sMesh.r[p], sMesh.z[s]);
            }
         }
      }
   }
}

void addSecondBoundaryConditions(ShiftArrays &shifts, SplittingMesh &sMesh, SLAE &slae, BoundaryConditions &conditions, FunctionProblem &funcs)
{
   Matrices mat{ };

   auto &M = mat.M;
   auto &H = mat.H;
   auto &conds = conditions.conds;
   auto &thetta = funcs.thetta;
   auto &r = sMesh.r;
   auto &z = sMesh.z;

   const int rSize = sMesh.r.size();
   const int zSize = sMesh.z.size();
   const int countConds = conds.size();

   for (int k = 0; k < countConds; k++)
   {
      int typeCond = conds[k][0];

      if (typeCond == 2)
      {
         int numFunc = conds[k][1];
         int p0 = conds[k][2];
         int p1 = conds[k][3];
         int s0 = conds[k][4];
         int s1 = conds[k][5];

         if (p0 == p1)
         {
            int p = shifts.Ir[p0];

            s0 = shifts.Iz[s0];
            s1 = shifts.Iz[s1];

            for (int s = s0; s < s1; s++)
            {
               double hz = z[s + 1] - z[s];

               vector<double> thettaValue = { thetta[numFunc](r[p], z[s]), thetta[numFunc](r[p], z[s + 1]) };

               int globalNum1 = rSize * s + p;
               int globalNum2 = rSize * (s + 1) + p;

               double elemB1 = 0;
               double elemB2 = 0;
               for (int j = 0; j < 2; j++)
               {
                  elemB1 += hz * M[0][j] * thettaValue[j];
                  elemB2 += hz * M[1][j] * thettaValue[j];
               }

               slae.b[globalNum1] += elemB1;
               slae.b[globalNum2] += elemB2;
            }
         }
         else
         {
            int s = shifts.Iz[s0];

            p0 = shifts.Ir[p0];
            p1 = shifts.Ir[p1];

            for (int p = p0; p < p1; p++)
            {
               double hr = r[p + 1] - r[p];

               vector<double> thettaValue = { thetta[numFunc](r[p], z[s]), thetta[numFunc](r[p + 1], z[s]) };

               int globalNum1 = rSize * s + p;
               int globalNum2 = rSize * s + p + 1;

               double elemB1 = 0;
               double elemB2 = 0;
               for (int j = 0; j < 2; j++)
               {
                  elemB1 += hr * (r[p] * M[0][j] + hr * H[0][j]) * thettaValue[j];
                  elemB2 += hr * (r[p] * M[1][j] + hr * H[1][j]) * thettaValue[j];
               }

               slae.b[globalNum1] += elemB1;
               slae.b[globalNum2] += elemB2;
            }
         }
      }
   }
}