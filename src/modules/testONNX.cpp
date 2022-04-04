#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"
#include "tools/selection.hpp"
#include "tools/jetCorrections.hpp"
#include "tools/bTagWeights.hpp"
#include "tools/leptonSF.hpp"
#include "tools/leptonCorrections.hpp"
#include "tools/triggerSF.hpp"
#include "tools/mcWeights.hpp"

#include <TFile.h>
#include <TGraphErrors.h>
#include <TTreeReader.h>
#include <TF1.h>
#include <TVector3.h>
#include <TMath.h>
#include <TStyle.h>
#include <TNtuple.h>
#include <iostream>
#include <fstream>
#include <chrono>
#include <numeric>
#include <onnxruntime_cxx_api.h>

Config const &cfg=Config::get();

std::ostream &operator<<(std::ostream &os, const std::vector<int64_t> &input)
{
   for (auto const &i: input) {
      os << i << " ";
   }
   return os;
}

template <typename T>
T vectorProduct(const std::vector<T>& v)
{
   return std::accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}

extern "C"
void run(){
   
   std::string instanceName{"DNNregression"};
   std::string modelFilepath{(cfg.DNN_Path+"_op13.onnx").Data()};

   Ort::Env env(OrtLoggingLevel::ORT_LOGGING_LEVEL_INFO,instanceName.c_str());
   Ort::SessionOptions sessionOptions;
   sessionOptions.SetIntraOpNumThreads(1);
   
   Ort::Session session(env, modelFilepath.c_str(), sessionOptions);
   
   Ort::AllocatorWithDefaultOptions allocator;
   
   size_t numInputNodes = session.GetInputCount();
   size_t numOutputNodes = session.GetOutputCount();

   std::cout << "Number of Input Nodes: " << numInputNodes << std::endl;
   std::cout << "Number of Output Nodes: " << numOutputNodes << std::endl;

   const char* inputName = session.GetInputName(0, allocator);
   std::cout << "Input Name: " << inputName << std::endl;

   Ort::TypeInfo inputTypeInfo = session.GetInputTypeInfo(0);
   auto inputTensorInfo = inputTypeInfo.GetTensorTypeAndShapeInfo();

   ONNXTensorElementDataType inputType = inputTensorInfo.GetElementType();
   std::cout << "Input Type: " << inputType << std::endl;

   std::vector<int64_t> inputDims = inputTensorInfo.GetShape();
   inputDims[0] = 1;
   std::cout << "Input Dimensions: " << inputDims << std::endl;

   const char* outputName = session.GetOutputName(0, allocator);
   std::cout << "Output Name: " << outputName << std::endl;

   Ort::TypeInfo outputTypeInfo = session.GetOutputTypeInfo(0);
   auto outputTensorInfo = outputTypeInfo.GetTensorTypeAndShapeInfo();

   ONNXTensorElementDataType outputType = outputTensorInfo.GetElementType();
   std::cout << "Output Type: " << outputType << std::endl;

   std::vector<int64_t> outputDims = outputTensorInfo.GetShape();
   outputDims[0] = 1;
   std::cout << "Output Dimensions: " << outputDims << std::endl;
   
   
   std::vector<const char*> inputNames{inputName};
   std::vector<const char*> outputNames{outputName};
   std::vector<Ort::Value> inputTensors;
   std::vector<Ort::Value> outputTensors;
   
   // ~size_t inputTensorSize = 20;
   size_t inputTensorSize = vectorProduct(inputDims);
   std::vector<float> inputTensorValues(inputTensorSize);

   // ~size_t outputTensorSize = 2;
   size_t outputTensorSize = vectorProduct(outputDims);
   std::vector<float> outputTensorValues(outputTensorSize);
   
   
   Ort::MemoryInfo memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);
   
   inputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, inputTensorValues.data(), inputTensorSize, inputDims.data(), inputDims.size()));

   outputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, outputTensorValues.data(), outputTensorSize, outputDims.data(), outputDims.size()));
   
   session.Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);
   
   for (int i=0; i<20; i++){
      std::cout<<inputTensorValues.at(i)<<std::endl;
   }
   std::cout<<outputTensorValues.at(0)<<std::endl;
   std::cout<<outputTensorValues.at(1)<<std::endl;
   
   std::cout<<"-------------------------------------------------"<<std::endl;
   
   inputTensorValues.at(0) = 5.;
   
   for (int i=0; i<20; i++){
      std::cout<<inputTensorValues.at(i)<<std::endl;
   }
   
   session.Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);
   
   std::cout<<outputTensorValues.at(0)<<std::endl;
   std::cout<<outputTensorValues.at(1)<<std::endl;
   
   
}

