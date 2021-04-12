//Script to plot DNN resolution studies

#include "Config.hpp"
#include "tools/hist.hpp"
#include "tools/physics.hpp"
#include "tools/io.hpp"
#include "tools/weighters.hpp"

#include <TFile.h>
#include <TF1.h>
#include <TF2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TColor.h>
#include "cppflow/ops.h"
#include "cppflow/model.h"

#include <chrono>


Config const &cfg=Config::get();

extern "C"
void run()
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;
    float total=0;
    
    // ~auto input = cppflow::fill({1, 47}, 2.0f);
    // ~auto input = cppflow::fill({1, 47}, tensor);
    cppflow::model model("/home/home4/institut_1b/dmeuser/top_analysis/pyKeras_ttbar/TMVA_pyKeras_reg/h5_to_pb/model");
    // ~auto output = model({{"serving_default_dense_1_input:0", input}},{"StatefulPartitionedCall:0"});
    
    for (int y=0; y<1000;y++){
        std::vector<float> input_vec(47);
        std::vector<int64_t> shape (2);
        for(int i=0;i<47;i++){
            input_vec[i]=1.0*i;
        }
        shape[0]=1;
        shape[1]=47;
        auto tensor = cppflow::tensor(input_vec, shape);
        auto t1 = high_resolution_clock::now();
        auto output = model({{"serving_default_dense_1_input:0", tensor}},{"StatefulPartitionedCall:0"});
        auto t2 = high_resolution_clock::now();
        
        auto values = output[0].get_data<float>();
        // ~auto values = output.data<float>();
        
        // ~for (auto v : values) {
            // ~std::cout << v << std::endl;
        // ~}
        duration<double, std::milli> ms_double = t2 - t1;
        total+=ms_double.count();
        std::cout << ms_double.count() << "ms\n";
    }
    std::cout<<"mean="<<total/1000<<std::endl;
}
