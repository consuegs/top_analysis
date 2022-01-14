#include "mcWeights.hpp"
#include "Config.hpp"
#include <iostream>

mcWeights::mcWeights(const Systematic::Systematic& systematic, const bool &isMadGraph, const bool &isPythiaOnly) :
systematic_(systematic),
isMadGraph_(isMadGraph),
isPythiaOnly_(isPythiaOnly)
{
   std::cout<<"--- Beginning preparation of mcWeights\n";
   
   // Check if syst. variation should be used
  const Systematic::Type type = systematic_.type();
   if(std::find(Systematic::mcWeightTypes.begin(), Systematic::mcWeightTypes.end(), type) != Systematic::mcWeightTypes.end()){
      useNominal = false;
      std::cout<<"Use systematic of type: "<<Systematic::convertType(type)<<"\n";
      if(systematic.variation() == Systematic::up) upVariation = true;
      else if(systematic.variation() == Systematic::down) upVariation = false;
      else{
         std::cerr<<"ERROR in constructor of mcWeights! Systematic variation is invalid: "
              <<Systematic::convertVariation(systematic.variation())<<"\n...break\n"<<std::endl;
         exit(205);
      }
   }
   else {
      useNominal = true;
      std::cout<<"Do not apply systematic variation, use nominal MC weight\n";
   }
   
   // Define correct binNr for me scales (different between powheg and madgraph)
   if(isMadGraph_) meScaleBins = {4,7,2,3,5,9};
   else meScaleBins = {2,3,4,7,5,9};
   
   //At this point lumi weight is not yet defined
   preparedLumiWeight = false;

}

// define current dataset and derive lumi weight
void mcWeights::prepareLumiWeight(const Datasubset &dss, const float &lumi){
   preparedLumiWeight = true;
   lumiWeight = dss.xsec/dss.getNgen_syst(systematic_)*lumi;
}


const float mcWeights::getMCweight(const float &nominalWeight, const std::vector<float> &w_pdf, const std::vector<float> &w_ps, const std::vector<float> &w_bFrag)
{
   if(useNominal) return nominalWeight;
   else {      // Compared to normalization in dataset.cpp Bins are shifted by -1 since entry in w_pdf[0] is in BinNR=1 for summed Hist
      switch(systematic_.type()){
         case Systematic::meFacScale:
            if (isPythiaOnly_) return nominalWeight;  //no PDF weights available for pythiaOnly
            else return w_pdf[(upVariation)? meScaleBins[0]-1 : meScaleBins[1]-1]*nominalWeight;
            break;
         case Systematic::meRenScale:
            if (isPythiaOnly_) return nominalWeight;  //no PDF weights available for pythiaOnly
            else return w_pdf[(upVariation)? meScaleBins[2]-1 : meScaleBins[3]-1]*nominalWeight;
            break;
         case Systematic::meScale:
            if (isPythiaOnly_) return nominalWeight;  //no PDF weights available for pythiaOnly
            else return w_pdf[(upVariation)? meScaleBins[4]-1 : meScaleBins[5]-1]*nominalWeight;
            break;
         case Systematic::alphasPdf:
            if (isPythiaOnly_) return nominalWeight;  //no PDF weights available for pythiaOnly
            else return w_pdf[(upVariation)? 112-1 : 111-1]*nominalWeight;
            break;
         case Systematic::psISRScale:
            return w_ps[(upVariation)? 27-1 : 26-1]*nominalWeight;
            break;
         case Systematic::psFSRScale:
            return w_ps[(upVariation)? 5-1 : 4-1]*nominalWeight;
            break;
         case Systematic::bFrag:
            return w_bFrag[(upVariation)? 1-1 : 3-1]*nominalWeight;
            break;
         case Systematic::bSemilep:
            return w_bFrag[(upVariation)? 5-1 : 6-1]*nominalWeight;
            break;
         default:
            std::cout<<"Error in mcWeights: Weight for "<<Systematic::convertType(systematic_.type())<<" not found"<<std::endl;
            exit(200);
            break;
      }
   }
}

const float mcWeights::getMCweight_lumiWeighted(const float &nominalWeight, const std::vector<float> &w_pdf, const std::vector<float> &w_ps, const std::vector<float> &w_bFrag)
{
   if (preparedLumiWeight){
      return this->getMCweight(nominalWeight,w_pdf,w_ps,w_bFrag)*lumiWeight;
   }
   else{
      std::cout<<"Error in mcWeights: lumiWeights is not defined!!"<<std::endl;
      exit(222);
   }
}
