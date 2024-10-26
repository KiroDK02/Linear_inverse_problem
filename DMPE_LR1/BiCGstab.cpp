#include "functions.h"

void BiCGStabWithLUPrecond(SLAE &slae, MatrixM &M, Vectors &vec, double eps, int maxIterations)
{
   const int N = slae.A.di.size();

   vec.r.resize(N);
   vec.r0.resize(N);
   vec.z.resize(N);
   vec.p.resize(N);
   vec.ay.resize(N);
   vec.az.resize(N);
   vec.s.resize(N);
   vec.y.resize(N);
   vec.x.resize(N);
   vec.b = slae.b;

   allocMemoryForM(slae.A, M);
   LUDecomp(slae, M);

   matrixVectorMult(slae.A, vec.x, vec.az);
   calculateResidual(vec.r0, vec.b, vec.az);
   double normB{ norm(vec.b) };
   double relativeResidual{ norm(vec.r0) / normB };
   vec.r = vec.r0;
   vec.p = vec.r0;
   double scalarRR0{ scalarProduct(vec.r, vec.r0) };

   for (int k{ 1 }; relativeResidual > eps && k < maxIterations; ++k)
   {
      calculateSLAE(slae, M, vec.y, vec.p);

      matrixVectorMult(slae.A, vec.y, vec.ay);
      double alpha{ scalarRR0 / scalarProduct(vec.ay, vec.r0) };

      additionVectorsWithCoef(vec.s, vec.r, vec.ay, 1, -alpha);

      calculateSLAE(slae, M, vec.z, vec.s);

      matrixVectorMult(slae.A, vec.z, vec.az);
      double omega{ scalarProduct(vec.az, vec.s) / scalarProduct(vec.az, vec.az) };

      additionVectorsWithCoef(vec.x, vec.x, vec.y, vec.z, 1, alpha, omega);

      additionVectorsWithCoef(vec.r, vec.s, vec.az, 1, -omega);

      double scalarRR0Prev{ scalarRR0 };
      scalarRR0 = scalarProduct(vec.r, vec.r0);
      double beta{ (alpha * scalarRR0) / (omega * scalarRR0Prev) };

      additionVectorsWithCoef(vec.p, vec.r, vec.p, vec.ay, 1, beta, -beta * omega);

      relativeResidual = norm(vec.r) / normB;

      std::cout << "Iteration #" << k << " , residual: " << relativeResidual << '\n';
   }

   matrixVectorMult(slae.A, vec.x, vec.az);
   calculateResidual(vec.r, vec.b, vec.az);
   relativeResidual = norm(vec.r) / normB;

   std::swap(slae.q, vec.x);
   std::cout << "Real residual: " << relativeResidual << '\n';
}

void allocMemoryForM(SparseMatrix &A, MatrixM &M)
{
   const int N = A.di.size();

   int dim{ A.ig[N] };

   M.di.resize(N);
   M.al.resize(dim);
   M.au.resize(dim);
}

void LUDecomp(SLAE &slae, MatrixM &M)
{
   const int N = slae.A.di.size();

   auto &ia = slae.A.ig;
   auto &ja = slae.A.jg;
   auto &au = slae.A.ggu;
   auto &al = slae.A.ggl;
   auto &di = slae.A.di;

   for (int i{ 0 }; i < N; ++i)
   {
      int i0{ ia[i] };
      int i1{ ia[i + 1] };

      int k{ i0 };
      double sumDi{ 0 };
      for (k; k < i1; k++)
      {
         int ki{ i0 };
         int j{ ja[k] };
         int j0{ ia[j] };
         int j1{ ia[j + 1] };
         int kj{ j0 };
         double sumAl{ 0 };
         double sumAu{ 0 };
         while (ki < k && kj < j1)
         {
            int ji{ ja[ki] };
            int jj{ ja[kj] };
            if (ji > jj) kj++;
            else if (ji < jj) ki++;
            else
            {
               sumAl += M.al[ki] * M.au[kj];
               sumAu += M.al[kj] * M.au[ki];
               ki++;
               kj++;
            }
         }

         M.al[k] = (al[k] - sumAl);
         M.au[k] = (au[k] - sumAu) / M.di[j];
         sumDi += M.al[k] * M.au[k];
      }

      M.di[i] = di[i] - sumDi;
   }
}

void calculateY(SLAE &slae, const MatrixM &M, std::vector<double> &y, const std::vector<double> &b)
{
   const int N = slae.A.di.size();

   auto &ia = slae.A.ig;
   auto &ja = slae.A.jg;
   auto &al = slae.A.ggl;

   for (int i{ 0 }; i < N; ++i)
   {
      int i0{ ia[i] };
      int i1{ ia[i + 1] };

      int k{ i0 };
      double s{ 0 };
      for (k; k < i1; ++k)
      {
         int j{ ja[k] };
         s += M.al[k] * y[j];
      }
      y[i] = (b[i] - s) / M.di[i];
   }
}


void calculateX(SLAE &slae, const MatrixM &M, std::vector<double> &x, const std::vector<double> &y)
{
   const int N = slae.A.di.size();
   
   auto &ia = slae.A.ig;
   auto &ja = slae.A.jg;
   auto &au = slae.A.ggu;
   
   x = y;

   for (int i{ N - 1 }; i >= 0; --i)
   {
      int i0{ ia[i] };
      int i1{ ia[i + 1] };

      int k{ i1 - 1 };
      double xi{ x[i] };
      for (k; k >= i0; --k)
      {
         int j{ ja[k] };
         x[j] -= M.au[k] * xi;
      }
      x[i] = xi;
   }
}

void calculateSLAE(SLAE &slae, const MatrixM &M, std::vector<double> &x, const std::vector<double> &y)
{
   calculateY(slae, M, x, y);
   calculateX(slae, M, x, x);
}

double scalarProduct(const std::vector<double> &x, const std::vector<double> &y)
{
   int N{ static_cast<int>(x.size()) };
   double sum{ 0 };

   for (int i{ 0 }; i < N; ++i)
   {
      sum += x[i] * y[i];
   }

   return sum;
}

void additionVectorsWithCoef(std::vector<double> &res, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, double coef1, double coef2, double coef3)
{
   int N{ static_cast<int>(x.size()) };

   for (int i{ 0 }; i < N; ++i)
   {
      res[i] = coef1 * x[i] + coef2 * y[i] + coef3 * z[i];
   }
}

void additionVectorsWithCoef(std::vector<double> &res, const std::vector<double> &x, const std::vector<double> &y, double coef1, double coef2)
{
   int N{ static_cast<int>(x.size()) };

   for (int i{ 0 }; i < N; ++i)
   {
      res[i] = coef1 * x[i] + coef2 * y[i];
   }
}

void calculateResidual(std::vector<double> &r, const std::vector<double> &b, const std::vector<double> &ax)
{
   int N{ static_cast<int>(r.size()) };

   for (int i{ 0 }; i < N; ++i)
   {
      r[i] = b[i] - ax[i];
   }
}

double norm(const std::vector<double> &x)
{
   int N{ static_cast<int>(x.size()) };
   double sum{ 0 };

   for (int i{ 0 }; i < N; ++i)
   {
      sum += x[i] * x[i];
   }

   return std::sqrt(sum);
}

void matrixVectorMult(SparseMatrix &A, const std::vector<double> &x, std::vector<double> &b)
{
   const int N = A.di.size();

   auto &ia = A.ig;
   auto &ja = A.jg;
   auto &al = A.ggl;
   auto &au = A.ggu;
   auto &di = A.di;

   for (int i{ 0 }; i < N; ++i)
   {
      int i0{ ia[i] };
      int i1{ ia[i + 1] };

      int k{ i0 };
      double s{ 0 };
      b[i] = di[i] * x[i];
      for (k; k < i1; ++k)
      {
         int j{ ja[k] };

         s += al[k] * x[j];
         b[j] += au[k] * x[i];
      }
      b[i] += s;
   }
}