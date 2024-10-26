#include "functions.h"

void inputMesh(Mesh &mesh)
{
   std::ifstream m("mesh.txt");

   int rMsize = 0;
   int zMsize = 0;
   int numElements = 0;

   m >> rMsize;
   mesh.rMain.resize(rMsize);

   for (int i = 0; i < rMsize; i++)
      m >> mesh.rMain[i];

   m >> zMsize;
   mesh.zMain.resize(zMsize);

   for (int i = 0; i < zMsize; i++)
      m >> mesh.zMain[i];

   m >> numElements;
   mesh.areasMesh.resize(numElements);

   for (int i = 0; i < numElements; i++)
   {
      mesh.areasMesh[i].resize(5);

      int numArea = 0;
      int r0 = 0;
      int r1 = 0;
      int z0 = 0;
      int z1 = 0;

      m >> numArea;
      m >> r0 >> r1;
      m >> z0 >> z1;

      mesh.areasMesh[i][0] = numArea - 1;
      mesh.areasMesh[i][1] = r0 - 1;
      mesh.areasMesh[i][2] = r1 - 1;
      mesh.areasMesh[i][3] = z0 - 1;
      mesh.areasMesh[i][4] = z1 - 1;
   }

   m.close();
}

// Считываем разбиение сетки (считывать будем для каждой оси просто передавая нужные массивы)
void inputSplittingMesh(vector<double> &coordAxisOfMesh, vector<double> &coordAxisOfSplitMesh, vector<int> &shiftArray, std::ifstream &split)
{
   const int AxisCoordSize = coordAxisOfMesh.size();

   shiftArray.resize(AxisCoordSize);

   int nk = 0;
   coordAxisOfSplitMesh.resize(1, coordAxisOfMesh[0]);

   for (int i = 0, j = 1; i < AxisCoordSize - 1; i++, j++)
   {
      int countIntervals = 0;

      double coef = 0;
      double step = 0;

      split >> countIntervals >> coef;
      nk += countIntervals;
      coordAxisOfSplitMesh.resize(nk + 1);

      if (coef != 1)
      {
         double sumProgression = (pow(coef, countIntervals) - 1.) / (coef - 1.);
         step = (coordAxisOfMesh[i + 1] - coordAxisOfMesh[i]) / sumProgression;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            coordAxisOfSplitMesh[j] = coordAxisOfMesh[i] + step * (pow(coef, jk) - 1.) / (coef - 1.);
      }
      else
      {
         step = (coordAxisOfMesh[i + 1] - coordAxisOfMesh[i]) / countIntervals;

         int jk = 1;
         for (j; j < nk; j++, jk++)
            coordAxisOfSplitMesh[j] = coordAxisOfMesh[i] + step * jk;
      }

      coordAxisOfSplitMesh[j] = coordAxisOfMesh[i + 1];
      shiftArray[i + 1] = j;
   }
}

void inputBoundaryCondtitions(BoundaryConditions &conditions)
{
   auto &conds = conditions.conds;
   int countConditions = 0;

   std::ifstream sConditions("conditions.txt");

   sConditions >> countConditions;
   conds.resize(countConditions);

   for (int i = 0; i < countConditions; i++)
   {
      conds[i].resize(6);

      int typeCond = 0;
      int numFunction = 0;

      int r0 = 0;
      int r1 = 0;
      int z0 = 0;
      int z1 = 0;

      sConditions >> typeCond >> numFunction;
      sConditions >> r0 >> r1;
      sConditions >> z0 >> z1;

      conds[i][0] = typeCond;
      conds[i][1] = numFunction - 1;
      conds[i][2] = r0 - 1;
      conds[i][3] = r1 - 1;
      conds[i][4] = z0 - 1;
      conds[i][5] = z1 - 1;
   }

   sConditions.close();
}