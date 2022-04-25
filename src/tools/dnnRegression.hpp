#ifndef DNNREGRESSION_HPP__
#define DNNREGRESSION_HPP__

#include <onnxruntime_cxx_api.h>

/// Class for application DNN regression
class DNNregression{

public:

   /// Constructor
   DNNregression(std::string modelFilepath);

   /// Destructor
   ~DNNregression(){}

   /// Get output of DNN regression
   void evaluate(const std::vector<float> &inputVec, std::vector<float> &outputVec);

private:
   
   Ort::Session* session;
   
   std::vector<const char*> inputNames;
   std::vector<const char*> outputNames;
   
   std::vector<float> inputTensorValues;
   std::vector<float> outputTensorValues;
   
   std::vector<Ort::Value> inputTensors;
   std::vector<Ort::Value> outputTensors;
   
};




#endif /* DNNREGRESSION_HPP__ */
