#include "OptimizerConstraints.h"

#include <math.h>

#define   CONSTRAIN_NONE    0
#define   CONSTRAIN_MIN     1
#define   CONSTRAIN_MAX     2
#define   CONSTRAIN_RANGE   3

void OptimizerInterface::Dimension(int parameters)
   {
   point.Dimension(parameters);
   min.Dimension(parameters);
   max.Dimension(parameters);

   constraints.Dimension(parameters);
   constraints.Zero();
   }

void OptimizerInterface::SetMin(int parameter, double value)
   {
   constraints[parameter] |= CONSTRAIN_MIN;
   min[parameter] = value;
   }

void OptimizerInterface::SetMax(int parameter, double value)
   {
   constraints[parameter] |= CONSTRAIN_MAX;
   max[parameter] = value;
   }

void OptimizerInterface::SetRange(int parameter, double MIN, double MAX)
   {
   constraints[parameter] = CONSTRAIN_RANGE;
   min[parameter] = MIN;
   max[parameter] = MAX;
   }

void OptimizerInterface::Fix(int parameter, double value)
   {
   constraints[parameter] = CONSTRAIN_RANGE;
   min[parameter] = max[parameter] = value;
   }

void OptimizerInterface::ClearConstraints(int parameter)
   {
   constraints[parameter] = CONSTRAIN_NONE;
   }

void OptimizerInterface::ClearConstraints()
   {
   constraints.Zero();
   }

void OptimizerInterface::SetObjectiveFunction(ObjectiveFunction & function)
   {
   f = &function;
   }

double OptimizerInterface::Evaluate(Vector & vector)
   {
   Translate(vector, point);

   return f->Evaluate(point);
   }

void OptimizerInterface::Translate(Vector & unconstrained, Vector & constrained)
   {
   constrained.Dimension(constraints.Length());

   for (int i = 0, j = 0; i < constraints.Length(); i++)
      switch (constraints[i])
         {
         case CONSTRAIN_NONE :
            constrained[i] = unconstrained[j++];
            break;
         case CONSTRAIN_MIN :
            constrained[i] = min[i] + exp(unconstrained[j++]);
            break;
         case CONSTRAIN_MAX :
            constrained[i] = max[i] - exp(unconstrained[j++]);
            break;
         case CONSTRAIN_RANGE :
            if (min[i] == max[i])
               constrained[i] = min[i];
            else
               {
               double x = unconstrained[j++];

               if (x >= 36)
                  constrained[i] = max[i];
               else
                  constrained[i] = min[i] + (max[i] - min[i]) * exp(x) / (1 + exp(x));
               }
         }
   }

void OptimizerInterface::BackTranslate(Vector & constrained, Vector & unconstrained)
   {
   unconstrained.Dimension(constraints.Length());

   int j = 0;
   for (int i = 0; i < constraints.Length(); i++)
      switch (constraints[i])
         {
         case CONSTRAIN_NONE :
            unconstrained[j++] = constrained[i];
            break;
         case CONSTRAIN_MIN :
            assert(constrained[i] >= min[i]);
            unconstrained[j++] = log(constrained[i] + min[i] + 1e-16);
            break;
         case CONSTRAIN_MAX :
            assert(constrained[i] <= max[i]);
            unconstrained[j++] = exp(max[i] - constrained[i] + 1e-16);
            break;
         case CONSTRAIN_RANGE :
            if (min[i] == max[i])
               assert(constrained[i] == min[i]);
            else
               {
               assert(constrained[i] >= min[i]);
               assert(constrained[i] <= max[i]);

               double x = (constrained[i] - min[i]) / (max[i] - min[i]);

               if (x >= 0.999999999) x = 0.999999999;
               if (x <= 1e-16)       x = 1e-16;

               unconstrained[j++] = log(x/(1-x));
               }
         }

   unconstrained.Dimension(j);
   }

int OptimizerInterface::CountParameters()
   {
   return constraints.Length();
   }

int OptimizerInterface::CountFreeParameters()
   {
   int parameters = constraints.Length();

   for (int i = 0; i < constraints.Length(); i++)
      if (constraints[i] == CONSTRAIN_RANGE && min[i] == max[i])
         parameters--;

   return parameters;
   }
