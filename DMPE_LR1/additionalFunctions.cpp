#include "functions.h"

double NumericalValueAtPoint(SplittingMesh &sMesh, vector<double> &q, double r, double z)
{
   auto &rCoords = sMesh.r;
   auto &zCoords = sMesh.z;

   const int rSize = rCoords.size();
   const int zSize = zCoords.size();

   int begR = 0;
   int endR = rSize - 1;
   int begZ = 0;
   int endZ = zSize - 1;

   while (!(rCoords[begR] <= r && r <= rCoords[begR + 1]))
   {
      int ind = (begR + endR) / 2;

      if (rCoords[ind] < r)
         begR = ind;
      else
         endR = ind;
   }

   while (!(zCoords[begZ] <= z && z <= zCoords[begZ + 1]))
   {
      int ind = (begZ + endZ) / 2;

      if (zCoords[ind] < z)
         begZ = ind;
      else
         endZ = ind;
   }

   double r0 = rCoords[begR];
   double r1 = rCoords[begR + 1];
   double z0 = zCoords[begZ];
   double z1 = zCoords[begZ + 1];

   double hr = r1 - r0;
   double hz = z1 - z0;

   double psi1 = (r1 - r) / hr * (z1 - z) / hz;
   double psi2 = (r - r0) / hr * (z1 - z) / hz;
   double psi3 = (r1 - r) / hr * (z - z0) / hz;
   double psi4 = (r - r0) / hr * (z - z0) / hz;

   vector<int> globalNums = { begZ * rSize + begR, begZ * rSize + begR + 1,
                              (begZ + 1) * rSize + begR, (begZ + 1) * rSize + begR + 1 };

   double numericalValue = q[globalNums[0]] * psi1 +
                           q[globalNums[1]] * psi2 + 
                           q[globalNums[2]] * psi3 + 
                           q[globalNums[3]] * psi4;

   return numericalValue;
}

int mu(int i)
{
   return i % 2;
}

int nu(int i)
{
   return i / 2;
}

void setSizeMatrix(int n, int m, vector<vector<double>> &matrix)
{
   matrix.resize(n);

   for (int i = 0; i < n; i++)
      matrix[i].resize(m);
}

void setStringMatrixInZero(SparseMatrix &matrix, int numString)
{
   int i0 = matrix.ig[numString];
   int i1 = matrix.ig[numString + 1];
   const int countNonZeroElems = matrix.ig[matrix.di.size()];

   // Зануляем всю строку i в нижнем треугольнике
   for (i0; i0 < i1; i0++)
      matrix.ggl[i0] = 0;

   int j0 = matrix.ig[numString + 1];

   // Проходим по столбцам, начиная с (i + 1)го 
   // до последнего ненулевого элемента в верхнем треугольнике
   // Зануляем элемент, у которого номер строки совпал с numString
   for (j0; j0 < countNonZeroElems; j0++)
      if (matrix.jg[j0] == numString)
         matrix.ggu[j0] = 0;

   matrix.di[numString] = 1;
}