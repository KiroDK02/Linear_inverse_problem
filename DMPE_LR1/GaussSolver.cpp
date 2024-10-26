#include "functions.h"

void MethodGaussSolver(SLAEFull &slae)
{
   auto &A = slae.A;
   auto &x = slae.x;
   auto &b = slae.b;

   const int N = A.size();

   for (int i = 0; i < N; i++)
   {
      int ni = 0;

      searchLeadingElement(A, i, ni);

      A[i].swap(A[ni]);
      std::swap(b[i], b[ni]);

      for (int k = i + 1; k < N; k++)
      {
         double coef = A[k][i] / A[i][i];

         for (int j = i; j < N; j++)
            A[k][j] -= coef * A[i][j];

         b[k] -= coef * b[i];
      }
   }

   for (int i = N - 1; i >= 0; i--)
   {
      double sum = 0;

      if (abs(A[i][i]) < 1e-15)
      {
         std::cout << "The matrix is degenerate.";
         exit(-1);
      }

      for (int j = N - 1; j > i; j--)
         sum += A[i][j] * x[j];

      x[i] = (b[i] - sum) / A[i][i];
   }
}

void searchLeadingElement(vector<vector<double>> &A, int i, int &ni)
{
   const int N = A.size();

   double lead = A[i][i];
   ni = i;

   for (int k = i + 1; k < N; k++)
      if (abs(A[k][i]) > abs(lead))
      {
         lead = A[k][i];
         ni = k;
      }
}