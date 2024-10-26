#include "functions.h"

const double alpha = 1e-13;
const vector<vector<double>> E = { { 1, 0, 0 },
                                   { 0, 1, 0 },
                                   { 0, 0, 1 } };

void solutionInverseProblem(SplittingMesh &sMesh, vector<double> &q, FunctionProblem &funcs, vector<double> &I0, vector<double> &result)
{
   auto &IReal = funcs.I;

   Sources sources;
   Receivers receivers;
   SLAEFull slae{ };

   const int countReceivers = receivers.M.size();
   const int m = sources.A.size();

   vector<double> IPrev = I0;
   vector<double> In = I0;
   vector<double> Inew(m);
   vector<double> UReal(countReceivers);
   vector<double> UIn(countReceivers);
   vector<double> w(countReceivers);
   
   setSizeMatrix(m, m, slae.A);
   slae.x.resize(m);
   slae.b.resize(m);
   result.resize(m);

   calcVoltageOnAllReceiver(UReal, sMesh, q, IReal, sources, receivers);
   calcVoltageOnAllReceiver(UIn, sMesh, q, I0, sources, receivers);

   setWeights(UReal, w);

   double functional = calcFunctional(UReal, UIn, w, In, IPrev);

   printf_s("I0 = (%.15lf, %.15lf, %.15lf)\nIReal = (%.15lf, %.15lf, %.15lf); functional = %.15lf\n", I0[0], I0[1], I0[2], IReal[0], IReal[1], IReal[2], functional);
   printf_s("n          I1                  I2                 I3             functional\n");

   for (int n = 0; n < 10000 && functional > 1e-14; n++)
   {
      buildSLAEForInverseProblem(UReal, slae, w, sMesh, q, In, IPrev, sources, receivers);
      MethodGaussSolver(slae);
      
      calcVecPlusCoefVec(Inew, In, 1, slae.x);

      std::swap(IPrev, In);
      std::swap(In, Inew);

      calcVoltageOnAllReceiver(UIn, sMesh, q, In, sources, receivers);

      functional = calcFunctional(UReal, UIn, w, In, IPrev);

      printf_s("%d; (%.15lf, %.15lf, %.15lf); %.15lf\n", n + 1, In[0], In[1], In[2], functional);
   }

   std::swap(In, result);
}

double calcFunctional(vector<double> &UReal, vector<double> UIn, vector<double> &w, vector<double> &In, vector<double> &IPrev)
{
   const int countReceivers = UReal.size();
   const int m = In.size();

   double functional = 0;

   for (int i = 0; i < countReceivers; i++)
   {
      double deltaU = (UIn[i] - UReal[i]);
      functional += w[i] * w[i] * deltaU * deltaU;
   }

   for (int i = 0; i < m; i++)
   {
      double InMinusIPrev = In[i] - IPrev[i];
      functional += alpha * InMinusIPrev * InMinusIPrev;
   }

   return functional;
}

void setWeights(vector<double> &UReal, vector<double> &w)
{
   const int countWeights = UReal.size();

   for (int i = 0; i < countWeights; i++)
      w[i] = 1. / UReal[i];
}

void buildSLAEForInverseProblem(vector<double> &Ureal, SLAEFull &slae, vector<double> &w, SplittingMesh &sMesh, vector<double> &q, vector<double> In, vector<double> IPrev, Sources &sources, Receivers &receivers)
{
   const int countReceivers = receivers.M.size();
   const int m = sources.A.size();

   double Up = 0;
   double Us = 0;
   double deltaU = 0;

   for (int p = 0; p < m; p++)
   {
      for (int s = 0; s < m; s++)
      {
         double sumAps = 0;

         for (int j = 0; j < countReceivers; j++)
         {
            Up = calcVoltageOnReceiverFromOneSource(sMesh, q, 1, sources, receivers, j, p);
            Us = calcVoltageOnReceiverFromOneSource(sMesh, q, 1, sources, receivers, j, s);

            sumAps += w[j] * w[j] * Up * Us;
         }

         slae.A[p][s] = sumAps + alpha * E[p][s];
      }

      double sumBp = 0;

      for (int j = 0; j < countReceivers; j++)
      {
         double UIn = calcVoltageOnReceiver(sMesh, q, In, sources, receivers, j);
         deltaU = UIn - Ureal[j];

         Up = calcVoltageOnReceiverFromOneSource(sMesh, q, 1, sources, receivers, j, p);

         sumBp -= w[j] * w[j] * deltaU * Up;
      }

      slae.b[p] = sumBp - alpha * (In[p] - IPrev[p]);
   }
}

void calcVoltageOnAllReceiver(vector<double> &U, SplittingMesh &sMesh, vector<double> &q, vector<double> &I, Sources &sources, Receivers &receivers)
{
   const int countReceivers = receivers.M.size();
   
   for (int i = 0; i < countReceivers; i++)
      U[i] = calcVoltageOnReceiver(sMesh, q, I, sources, receivers, i);
}

double calcVoltageOnReceiver(SplittingMesh &sMesh, vector<double> &q, vector<double> &I, Sources &sources, Receivers &receivers, int numReceiver)
{
   const int countSources = sources.A.size();

   double valueVoltage = 0;

   for (int i = 0; i < countSources; i++)
      valueVoltage += calcVoltageOnReceiverFromOneSource(sMesh, q, I[i], sources, receivers, numReceiver, i);

   return valueVoltage;
}

double calcVoltageOnReceiverFromOneSource(SplittingMesh &sMesh, vector<double> &q, double I, Sources &sources, Receivers &receivers, int numReceiver, int numSource)
{
   double VBiM = calcPotentialAtPoint(sMesh, q, receivers.M[numReceiver], sources.B[numSource]);
   double VAiM = calcPotentialAtPoint(sMesh, q, receivers.M[numReceiver], sources.A[numSource]);
   double VBiN = calcPotentialAtPoint(sMesh, q, receivers.N[numReceiver], sources.B[numSource]);
   double VAiN = calcPotentialAtPoint(sMesh, q, receivers.N[numReceiver], sources.A[numSource]);

//   double valueVoltage = I / (2 * PI) * ((1. / rBM - 1. / rAM) - (1. / rBN - 1. / rAN));
   double valueVoltage = I * ((VBiM + VAiM) - (VBiN + VAiN));

   return valueVoltage;
}

double calcPotentialAtPoint(SplittingMesh &sMesh, vector<double> &q, vector<double> &pointReceiver, vector<double> &pointSource)
{
   vector<double> vec{ pointReceiver[0] - pointSource[0],
                       pointReceiver[1] - pointSource[1], 
                       pointReceiver[2] - pointSource[2] };

   double r = calcEuclideanNorm(vec);

   double valuePotential = NumericalValueAtPoint(sMesh, q, r, 0);

   return valuePotential;
}