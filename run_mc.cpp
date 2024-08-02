#include "./mc_subroutine/mc_read_load_compute.hpp"
#include "./potentialFunction/potentialFunctionPrototype.hpp"

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cout << "wrong arguments" << std::endl;
        std::exit(2);
    }

    auto mcObj=mc_computation(std::string(argv[1]));
//
//    auto funcPtr= createPotentialFunction("V_inv_12_6","25,80,15,67,2");
//    funcPtr->init();
//    int N=2;
//    std::shared_ptr<double[]>d0Vec= std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
//    std::shared_ptr<double[]>d1Vec= std::shared_ptr<double[]>(new double[N-1], std::default_delete<double[]>());
//    d0Vec[0]=1;
//    d0Vec[1]=1;
//    d1Vec[0]=1;
//    double L=4;
//    double  val=(*funcPtr)(L,d0Vec.get(),d1Vec.get());
//    std::cout<<"val="<<val<<std::endl;





}