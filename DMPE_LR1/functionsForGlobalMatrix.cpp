#include "functions.h"

void generatePortraitOfSparseMatrix(SLAE &slae, int rSize, int zSize)
{
   auto &ig = slae.A.ig;
   auto &jg = slae.A.jg;

   const int N = rSize * zSize;

   vector<std::set<int>> list{ };
   
   list.resize(N);
   slae.A.di.resize(N);
   slae.b.resize(N);
   slae.q.resize(N);

   for (int s = 0; s < zSize - 1; s++)
      for (int p = 0; p < rSize - 1; p++)
      {
         vector <int> globalNum{ (s + 1) * rSize + p + 1, (s + 1) * rSize + p,
                                  s * rSize + p + 1, s * rSize + p };

         for (int i = 0; i < 4; i++)
         {
            int ind1 = globalNum[i];

            for (int j = i + 1; j < 4; j++)
            {
               int ind2 = globalNum[j];
               list[ind1].insert(ind2);
            }
         }
      }

   ig.resize(N + 1);

   ig[0] = 0;
   ig[1] = 0;

   for (int i = 0; i < N; i++)
      ig[i + 1] = ig[i] + list[i].size();

   int countNonZeroElem = ig[N];

   jg.resize(countNonZeroElem);
   slae.A.ggl.resize(countNonZeroElem);
   slae.A.ggu.resize(countNonZeroElem);

   for (int i = 0, k = 0; i < N; i++)
      for (auto j : list[i])
      {
         jg[k] = j;
         k++;
      }
}

void calcLocalStiffnessMatrix(double r0, double hr, double hz, vector<vector<double>> &stiffness)
{
   Matrices mat{ };

   auto &G = mat.G;
   auto &M = mat.M;
   auto &H = mat.H;

   double coef = (r0 + hr / 2.) / hr;

   for (int i = 0; i < 4; i++)
   {
      int mui = mu(i);
      int nui = nu(i);
      for (int j = 0; j < 4; j++)
      {
         int muj = mu(j);
         int nuj = nu(j);
         stiffness[i][j] = coef * G[mui][muj] * hz * M[nui][nuj] + G[nui][nuj] / hz * hr * (r0 * M[mui][muj] + hr * H[mui][muj]);
      }
   }
}

void calcLocalRightPart(vector<double> &valueF, double r0, double hr, double hz, vector<double> &result)
{
   Matrices mat{ };

   auto &M = mat.M;
   auto &H = mat.H;

   for (int i = 0; i < 4; i++)
   {
      int mui = mu(i);
      int nui = nu(i);
      double elemb = 0;

      for (int j = 0; j < 4; j++)
      {
         int muj = mu(j);
         int nuj = nu(j);

         elemb += hr * (r0 * M[mui][muj] + hr * H[mui][muj]) * hz * M[nui][nuj] * valueF[j];
      }

      result[i] = elemb;
   }
}

void addLocalMatrixToGlobalMatrix(vector<int> &globalNumbers, double sigma, vector<vector<double>> &localMatrix, SLAE &slae)
{
   for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
      {
         double elem = sigma * localMatrix[i][j];

         addLocalElemToGlobalMatrix(slae, elem, globalNumbers[i], globalNumbers[j]);
      }
}

void addLocalRightPartToGlobal(vector<int> &globalNumbers, vector<double> &localRightPart, SLAE &slae)
{
   auto &b = slae.b;

   for (int i = 0; i < 4; i++)
      b[globalNumbers[i]] += localRightPart[i];
}

void addLocalElemToGlobalMatrix(SLAE &slae, double elem, int i, int j)
{
   auto &A = slae.A;

   if (i == j)
      A.di[i] += elem;
   else
   {
      if (i > j)
      {
         int beg = A.ig[i];
         int end = A.ig[i + 1] - 1;

         while (A.jg[beg] != j)
         {
            int ind = (beg + end) / 2;

            if (A.jg[ind] < j)
               beg = ind + 1;
            else
               end = ind;
         }

         A.ggl[beg] += elem;
      }
      else
      {
         int beg = A.ig[j];
         int end = A.ig[j + 1] - 1;

         while (A.jg[beg] != i)
         {
            int ind = (beg + end) / 2;

            if (A.jg[ind] < i)
               beg = ind + 1;
            else
               end = ind;
         }

         A.ggu[beg] += elem;
      }
   }
}