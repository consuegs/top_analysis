#include "dnnRegression.hpp"
#include <iostream>
#include <numeric>

template <typename T>
T vectorProduct(const std::vector<T>& v)
{
   return std::accumulate(v.begin(), v.end(), 1, std::multiplies<T>());
}

DNNregression::DNNregression(std::string modelFilepath)
{
   modelFilepath = modelFilepath+".onnx";
   std::cout<<"--- Applying DNN regression based on "<<modelFilepath<<std::endl;
   
   std::string instanceName{"DNNregression"};

   Ort::Env env(OrtLoggingLevel::ORT_LOGGING_LEVEL_WARNING,instanceName.c_str());
   Ort::SessionOptions sessionOptions;
   sessionOptions.SetIntraOpNumThreads(1);

   session = new Ort::Session(env, modelFilepath.c_str(), sessionOptions);
   
   Ort::AllocatorWithDefaultOptions allocator;
   
   const char* inputName = session->GetInputName(0, allocator);
   const char* outputName = session->GetOutputName(0, allocator);
   
   Ort::TypeInfo inputTypeInfo = session->GetInputTypeInfo(0);
   Ort::TypeInfo outputTypeInfo = session->GetOutputTypeInfo(0);
   auto inputTensorInfo = inputTypeInfo.GetTensorTypeAndShapeInfo();
   auto outputTensorInfo = outputTypeInfo.GetTensorTypeAndShapeInfo();
   
   std::vector<int64_t> inputDims = inputTensorInfo.GetShape();
   std::vector<int64_t> outputDims = outputTensorInfo.GetShape();
   inputDims[0] = 1;
   outputDims[0] = 1;
   
   size_t inputTensorSize = vectorProduct(inputDims);
   size_t outputTensorSize = vectorProduct(outputDims);
   inputTensorValues.resize(inputTensorSize);
   outputTensorValues.resize(outputTensorSize);
   
   inputNames = {inputName};
   outputNames = {outputName};
   
   Ort::MemoryInfo memoryInfo = Ort::MemoryInfo::CreateCpu(OrtAllocatorType::OrtArenaAllocator, OrtMemType::OrtMemTypeDefault);

   inputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, inputTensorValues.data(), inputTensorSize, inputDims.data(), inputDims.size()));
   outputTensors.push_back(Ort::Value::CreateTensor<float>(memoryInfo, outputTensorValues.data(), outputTensorSize, outputDims.data(), outputDims.size()));

}

void DNNregression::evaluate(const std::vector<float> &inputVec, std::vector<float> &outputVec)
{
   inputTensorValues = inputVec;
   
   session->Run(Ort::RunOptions{nullptr}, inputNames.data(), inputTensors.data(), 1, outputNames.data(), outputTensors.data(), 1);
   outputVec = outputTensorValues;
}


