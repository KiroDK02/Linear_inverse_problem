#include "functions.h"

int numberOfCurrentArea(ShiftArrays &shifts, vector<vector<int>> &areasMesh, int p, int s, int &l)
{
   auto &Ir = shifts.Ir;
   auto &Iz = shifts.Iz;

   const int countAreas = areasMesh.size();
   const int L1 = l;

   for (l; l < countAreas; l++)
   {
      const int numArea = areasMesh[l][0];
      const int p0 = Ir[areasMesh[l][1]];
      const int p1 = Ir[areasMesh[l][2]];
      const int s0 = Iz[areasMesh[l][3]];
      const int s1 = Iz[areasMesh[l][4]];

      if (p0 <= p && p <= p1 && p0 <= (p + 1) && (p + 1) <= p1 &&
          s0 <= s && s <= s1 && s0 <= (s + 1) && (s + 1) <= s1)
         return numArea;
   }

   for (l = 0; l < L1; l++)
   {
      const int numArea = areasMesh[l][0];
      const int p0 = Ir[areasMesh[l][1]];
      const int p1 = Ir[areasMesh[l][2]];
      const int s0 = Iz[areasMesh[l][3]];
      const int s1 = Iz[areasMesh[l][4]];

      if (p0 <= p && p <= p1 && p0 <= (p + 1) && (p + 1) <= p1 &&
         s0 <= s && s <= s1 && s0 <= (s + 1) && (s + 1) <= s1)
         return numArea;
   }

   return -1;
}

bool isFictitiousNode(ShiftArrays &shifts, vector<vector<int>> &areasMesh, int p, int s, int &l)
{
   auto &Ir = shifts.Ir;
   auto &Iz = shifts.Iz;

   const int countAreas = areasMesh.size();
   const int L1 = l;

   for (l; l < countAreas; l++)
   {
      const int p0 = Ir[areasMesh[l][1]];
      const int p1 = Ir[areasMesh[l][2]];
      const int s0 = Iz[areasMesh[l][3]];
      const int s1 = Iz[areasMesh[l][4]];

      if (p0 <= p && p <= p1 &&
          s0 <= s && s <= s1)
         return false;
   }

   for (l = 0; l < L1; l++)
   {
      const int p0 = Ir[areasMesh[l][1]];
      const int p1 = Ir[areasMesh[l][2]];
      const int s0 = Iz[areasMesh[l][3]];
      const int s1 = Iz[areasMesh[l][4]];

      if (p0 <= p && p <= p1 &&
         s0 <= s && s <= s1)
         return false;
   }

   return true;
}