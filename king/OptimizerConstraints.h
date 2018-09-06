#ifndef __OPTIMIZER_INTERFACE_H__
#define __OPTIMIZER_INTERFACE_H__

#include "MathVector.h"
#include "IntArray.h"

class ObjectiveFunction
   {
   public:
      virtual ~ObjectiveFunction() { };

      virtual double Evaluate(Vector & v) = 0;
   };

class OptimizerInterface : public VectorFunc
   {
   public:
      virtual  double Evaluate(Vector & v);

      void     Dimension(int parameters);
      int      CountFreeParameters();
      int      CountParameters();

      void     ClearConstraints();

      void     SetMin(int parameter, double min);
      void     SetMax(int parameter, double max);
      void     SetRange(int parameter, double min, double max);
      void     Fix(int parameter, double value);
      void     ClearConstraints(int parameter);

      void     SetObjectiveFunction(ObjectiveFunction & f);

      void     Translate(Vector & unconstrained, Vector & constrained);
      void     BackTranslate(Vector & constrained, Vector & unconstrained);

   private:
      IntArray constraints;
      Vector   min, max, point;

      ObjectiveFunction * f;
   };


#endif


