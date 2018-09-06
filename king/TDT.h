#ifndef __TDT_h__
#define __TDT_h__

class TDT:public AssociationAnalysis{
   public:
      TDT(Pedigree & pedigree);
      ~TDT();
      IntArray MT, MNT, pooT[2], pooNT[2];
      void SetupGlobals();
      void pre_genome();
      void post_genome();
      void PrintScores();

      void Analyze();
      void AnalyzeX();

      void TDT_Preparation();
      void TDT_Analysis();
      void GEE_Analysis();
      void WGDT_Analysis();
      void TDT1P_Analysis();
      void GDT_PO_Analysis();
      void GDT_CP_Analysis();
      void GDT_Missing_Analysis();
      void PDT_Analysis();
      void EPDT_Analysis();
      void RDT_Analysis();
      void GDT_Analysis();
      void GDT_hetero_Analysis();
      void FCAT_Analysis();

      void GDT_AnalysisX();
      void GDT_PO_AnalysisX();

      int CheckTwin();
      IntArray TwinFlag_Fam;
      IntArray TwinFlag_ID;
      void WrongQLS_Analysis();
      void MQLS_Analysis();

      char relationship(int i, int j);
//      String relationSet;

      // unreleased
      void QLS_Analysis();
      void bQLS_Analysis();
      void tQLS_Analysis();
      void cQLS_Analysis();
      void QLS_hetero_Analysis();
      void LogisticScore_Analysis();
      void GDT_FO_Analysis();
      void GDT_MO_Analysis();

      void rareGDT_Analysis();
      void rareEDA_Analysis();
//      void AutosomalCheck();
};

#endif
