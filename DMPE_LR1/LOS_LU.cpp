#include "functions.h"

void LOS_LU(SLAE &slae, int maxIter, double eps)
{
   const int sizeMatrix = slae.A.di.size();
   double normb = 0;

   SLAE LU{ };

   vector<double> r1(sizeMatrix);
   vector<double> z1(sizeMatrix);
   vector<double> p1(sizeMatrix);
   vector<double> mult(sizeMatrix);
   vector<double> rk(sizeMatrix);
   vector<double> Ar(sizeMatrix);
   vector<double> p(sizeMatrix);

   calcLU(slae, LU);
   calcDiscrepancy(slae, r1, normb);

   calcY(LU, r1, r1);
   calcX(LU, r1, z1);

   multSparseMatrixToVector(slae.A, z1, p1);
   calcY(LU, p1, p1);

   double scalarr = scalarMultiply(r1, r1);
   double discrepancy = sqrt(scalarr / normb);

   int k = 1;

   for (k; k < maxIter && discrepancy > eps; k++)
   {
      std::cout << discrepancy << "\n";
      double scalarp = scalarMultiply(p1, p1);
      double alpha = scalarMultiply(p1, r1) / scalarp;

      calcVectorMultToCoef(z1, alpha, mult);
      calcVectorPlusVector(slae.q, mult, slae.q);

      calcVectorMultToCoef(p1, -alpha, mult);
      calcVectorPlusVector(r1, mult, r1);

      calcX(LU, r1, rk);
      multSparseMatrixToVector(slae.A, rk, Ar);
      calcY(LU, Ar, p);

      double betta = -scalarMultiply(p1, p) / scalarp;

      calcVectorMultToCoef(z1, betta, mult);
      calcVectorPlusVector(rk, mult, z1);

      calcVectorMultToCoef(p1, betta, mult);
      calcVectorPlusVector(p, mult, p1);

      discrepancy = sqrt(scalarMultiply(r1, r1) / normb);
   }

   normb = 0;

   calcDiscrepancy(slae, r1, normb);
   
   discrepancy = sqrt(scalarMultiply(r1, r1) / normb);

   std::cout << discrepancy << "\n";
}

void calcDiscrepancy(SLAE &slae, vector<double> &discrepancy, double &normb)
{
   auto &matrix = slae.A;
   auto &b = slae.b;
   auto &q = slae.q;

   const int sizeMatrix = matrix.di.size();

   for (int i = 0; i < sizeMatrix; i++)
   {
      normb += b[i] * b[i];
      discrepancy[i] = b[i] - matrix.di[i] * q[i];

      int i0 = matrix.ig[i];
      int i1 = matrix.ig[i + 1];

      for (i0; i0 < i1; i0++)
      {
         int j = matrix.jg[i0];

         discrepancy[i] -= matrix.ggl[i0] * q[j];
         discrepancy[j] -= matrix.ggu[i0] * q[i];
      }
   }
}