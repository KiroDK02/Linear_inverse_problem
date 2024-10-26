#include "structures.h"

struct FunctionProblem
{
   vector<double> I
   {
      5, 6, 8
   };

   vector<double> sigma
   {
      1
   };

   vector<std::function<double(double, double)>> f
   {
      [](double r, double z) { return 0; }
   };

   vector<std::function<double(double, double)>> ug
   {
      [](double r, double z) { return 0; },
      [](double r, double z) { return 0; },
      [](double r, double z) { return 0; },
      [](double r, double z) { return z; }
   };

   vector<std::function<double(double, double)>> thetta
   {
      [](double r, double z) { return 0; },
      [](double r, double z) { return 0; }
   };
};

// Решение обратной задачи
void solutionInverseProblem(SplittingMesh &sMesh, vector<double> &q, FunctionProblem &funcs, vector<double> &I0, vector<double> &result);
double calcFunctional(vector<double> &UReal, vector<double> UIn, vector<double> &w, vector<double> &In, vector<double> &IPrev);
void setWeights(vector<double> &UReal, vector<double> &w);
void buildSLAEForInverseProblem(vector<double> &Ureal, SLAEFull &slae, vector<double> &w, SplittingMesh &sMesh, vector<double> &q, vector<double> In, vector<double> IPrev, Sources &sources, Receivers &receivers);
void calcVoltageOnAllReceiver(vector<double> &U, SplittingMesh &sMesh, vector<double> &q, vector<double> &I, Sources &sources, Receivers &receivers);
double calcVoltageOnReceiver(SplittingMesh &sMesh, vector<double> &q, vector<double> &I, Sources &sources, Receivers &receivers, int numReceiver);
double calcVoltageOnReceiverFromOneSource(SplittingMesh &sMesh, vector<double> &q, double I, Sources &sources, Receivers &receivers, int numReceiver, int numSource);
double calcPotentialAtPoint(SplittingMesh &sMesh, vector<double> &q, vector<double> &pointReceiver, vector<double> &pointSource);

// Функции считывания
void inputMesh(Mesh &mesh);
void inputSplittingMesh(vector<double> &coordAxisOfMesh, vector<double> &coordAxisOfSplitMesh, vector<int> &shiftArray, std::ifstream &split);
void inputBoundaryCondtitions(BoundaryConditions &conditions);

// Функции для работы с расчетной областью
int numberOfCurrentArea(ShiftArrays &shifts, vector<vector<int>> &areasMesh, int p, int s, int &l);
bool isFictitiousNode(ShiftArrays &shifts, vector<vector<int>> &areasMesh, int p, int s, int &l);

// Функции сборки глобальной матрицы и вектора правой части
void calcGlobalMatrix(SLAE &slae, SplittingMesh &sMesh, ShiftArrays &shifts, Mesh &mesh, FunctionProblem &funcs, BoundaryConditions &conditions);
void addFirstBoundaryConditions(ShiftArrays &shifts, SplittingMesh &sMesh, SLAE &slae, BoundaryConditions &conditions, FunctionProblem &funcs);
void addSecondBoundaryConditions(ShiftArrays &shifts, SplittingMesh &sMesh, SLAE &slae, BoundaryConditions &conditions, FunctionProblem &funcs);

// Функции для глобальной матрицы и вектора правой части
void generatePortraitOfSparseMatrix(SLAE &slae, int rSize, int zSize);
void calcLocalStiffnessMatrix(double r0, double hr, double hz, vector<vector<double>> &stiffness);
void calcLocalRightPart(vector<double> &valueF, double r0, double hr, double hz, vector<double> &result);
void addLocalMatrixToGlobalMatrix(vector<int> &globalNumbers, double sigma, vector<vector<double>> &localMatrix, SLAE &slae);
void addLocalElemToGlobalMatrix(SLAE &slae, double elem, int i, int j);
void addLocalRightPartToGlobal(vector<int> &globalNumbers, vector<double> &localRightPart, SLAE &slae);

// Предобуславливание матрицы
void calcLU(SLAE &slae, SLAE &LU);
void calcY(SLAE &LU, vector<double> &b, vector<double> &y);
void calcX(SLAE &LU, vector<double> &y, vector<double> &x);

// Решатель LOS для СЛАУ
void LOS_LU(SLAE &slae, int maxIter, double eps);
void calcDiscrepancy(SLAE &slae, vector<double> &discrepancy, double &normb);
void initFirstColumn(vector<vector<double>> &V, vector<double> &r);
void initSizeOfMatrix(vector<vector<double>> &matrix, int rows, int columns);
void initRotationMatrix(vector<vector<double>> &R, vector<vector<double>> &H, vector<vector<double>> &H0, int numColumn);
void doUpperTriangle(vector<vector<double>> &H, vector<double> &d, int m);
void calcColumnH(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &w, int mu);
void calcColumnV(vector<vector<double>> &V, double elemH, vector<double> &vAdditional, int mu);
void getColumnOfMatrix(vector<vector<double>> &matrix, vector<double> &column, int numColumn);
void calcAdditionalVector(vector<vector<double>> &H, vector<vector<double>> &V, vector<double> &vAdditional, int mu);

// Решатель GMRES для СЛАУ
int GMRES(SLAE &slae, int m, int maxIter, double eps);
void calcVectorDiscrepancy(SLAE &slae, vector<double> &discrepancy);

void BiCGStabWithLUPrecond(SLAE &slae, MatrixM &M, Vectors &vec, double eps, int maxIterations);
void matrixVectorMult(SparseMatrix &A, const std::vector<double> &x, std::vector<double> &b);
void allocMemoryForM(SparseMatrix &A, MatrixM &M);
void LUDecomp(SLAE &slae, MatrixM &M);
void calculateY(SLAE &slae, const MatrixM &M, std::vector<double> &y, const std::vector<double> &b);
void calculateX(SLAE &slae, const MatrixM &M, std::vector<double> &x, const std::vector<double> &y);
void calculateSLAE(SLAE &slae, const MatrixM &M, std::vector<double> &x, const std::vector<double> &y);
void inputVectors(std::vector<double> &x, std::vector<double> &b);
double scalarProduct(const std::vector<double> &x, const std::vector<double> &y);
void additionVectorsWithCoef(std::vector<double> &res, const std::vector<double> &x, const std::vector<double> &y, double coef1, double coef2);
void additionVectorsWithCoef(std::vector<double> &res, const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &z, double coef1, double coef2, double coef3);
void calculateResidual(std::vector<double> &r, const std::vector<double> &b, const std::vector<double> &ax);
double norm(const std::vector<double> &x);

// Решатель метод Гаусса для СЛАУ
void MethodGaussSolver(SLAEFull &slae);
void searchLeadingElement(vector<vector<double>> &A, int i, int &ni);

// Алгоритмы линейной алгебры
void multSparseMatrixToVector(SparseMatrix &A, vector<double> &vec, vector<double> &result);
void multMatrixToVector(vector<vector<double>> &matrix, vector<double> &vec, vector<double> &result);
void multLowTriangleToVector(SparseMatrix &matrix, vector<double> &vec, vector<double> &result);
double scalarMultiply(vector<double> &vector1, vector<double> &vector2);
void calcVectorMultToCoef(vector <double> &vec, double coef, vector<double> &result);
void calcVectorPlusVector(vector<double> &vector1, vector<double> &vector2, vector<double> &result);
double calcEuclideanNorm(vector<double> &vec);
void calcVecPlusCoefVec(vector<double> &result, vector<double> &vec1, double coef, vector<double> &vec2);
void calcSlaeWithUpTriangle(vector<vector<double>> &H, vector<double> &x, vector<double> &b, int sizeSlae);
void clearMatrix(vector<vector<double>> &matrix);
void equalMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2);
void multMatrixToMatrix(vector<vector<double>> &matrix1, vector<vector<double>> &matrix2, vector<vector<double>> &matrixResult);

// Доп функции
double NumericalValueAtPoint(SplittingMesh &sMesh, vector<double> &q, double r, double z);
int mu(int i);
int nu(int i);
void setSizeMatrix(int n, int m, vector<vector<double>> &matrix);
void setStringMatrixInZero(SparseMatrix &matrix, int numString);