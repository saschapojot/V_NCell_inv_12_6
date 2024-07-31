//
// Created by polya on 7/19/24.
//

#include "mc_read_load_compute.hpp"





void mc_computation::execute_mc(const double& L,const std::shared_ptr<double[]>& d0Vec, const std::shared_ptr<double[]>& d1Vec, const size_t & loopInit, const size_t & flushNum){

    double LCurr = L;
    std::shared_ptr<double[]> d0VecCurr=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

    std::shared_ptr<double[]> d1VecCurr=std::shared_ptr<double[]>(new double[N-1], std::default_delete<double[]>());

    std::memcpy(d0VecCurr.get(),d0Vec.get(),N*sizeof (double ));
    std::memcpy(d1VecCurr.get(),d1Vec.get(),(N-1)*sizeof (double ));

    //initialize next values
    double LNext;
    std::shared_ptr<double[]> d0VecNext=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
    std::shared_ptr<double[]> d1VecNext=std::shared_ptr<double[]>(new double[N-1], std::default_delete<double[]>());

    double UCurr;// = (*potFuncPtr)(LCurr, y0Curr, z0Curr, y1Curr);
    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    size_t loopStart = loopInit;
    for (size_t fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (size_t j = 0; j < loopToWrite; j++) {
            //propose a move
//            double LNext;
//            double y0Next;
//            double z0Next;
//            double y1Next;
//            double LReset;

            this->proposal(LCurr,d0VecCurr,d1VecCurr,LNext,d0VecNext,d1VecNext);
            double UNext;
            UCurr=(*potFuncPtr)(LCurr,d0VecCurr.get(),d1VecCurr.get());
            double r= acceptanceRatio(LCurr,d0VecCurr,d1VecCurr,UCurr,
                                      LNext,d0VecNext,d1VecNext,UNext);
            double u = distUnif01(e2);
            if (u <= r) {
                LCurr = LNext;
                std::memcpy(d0VecCurr.get(),d0VecNext.get(),N*sizeof (double ));
                std::memcpy(d1VecCurr.get(),d1VecNext.get(),(N-1)*sizeof (double ));


                UCurr = UNext;

            }//end of accept-reject
            U_dist_ptr[varNum*j+0]=UCurr;
            U_dist_ptr[varNum*j+1]=LCurr;
            for(int n=2;n<=2+N;n++){
                U_dist_ptr[varNum*j+n]=d0VecCurr[n-2];
            }
            for(int n=N+3;n<=2*N+1;n++){
                U_dist_ptr[varNum*j+n]=d1VecCurr[n-N-3];
            }

        }//end for loop
        size_t loopEnd = loopStart + loopToWrite - 1;
        std::string fileNameMiddle = "loopStart" + std::to_string(loopStart) + "loopEnd" + std::to_string(loopEnd);
        std::string out_U_distPickleFileName = this->U_dist_dataDir + "/" + fileNameMiddle + ".U_dist.csv";

        //save U_dist_ptr
        saveArrayToCSV(U_dist_ptr,varNum * loopToWrite,out_U_distPickleFileName,varNum);
        const auto tMCEnd{std::chrono::steady_clock::now()};
        const std::chrono::duration<double> elapsed_secondsAll{tMCEnd - tMCStart};
        std::cout << "loop " + std::to_string(loopStart) + " to loop " + std::to_string(loopEnd) + ": "
                  << elapsed_secondsAll.count() / 3600.0 << " h" << std::endl;

        loopStart = loopEnd + 1;
    }//end flush for loop

    std::cout<<"mc executed for "<<flushNum<<" flushes."<<std::endl;


}







void mc_computation::saveArrayToCSV(const std::shared_ptr<double[]>& array, const  int& arraySize, const std::string& filename, const int& numbersPerRow) {

    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    outFile<<"U,"<<"L";
    for(int i=0;i<N;i++){
        outFile<<",d0"+std::to_string(i);
    }
    for(int i=0;i<N-1;i++){
        outFile<<",d1"+std::to_string(i);
    }
    outFile<<"\n";
    for (int i = 0; i < arraySize; ++i) {
        outFile << array[i];
        if ((i + 1) % numbersPerRow == 0) {
            outFile << '\n';
        } else {
            outFile << ',';
        }
    }

    // If the last row isn't complete, end the line
    if (arraySize % numbersPerRow != 0) {
        outFile << '\n';
    }

    outFile.close();


}

void mc_computation::init_and_run(){
    this->execute_mc(LInit,d0VecInit,d1VecInit,loopLastFile+1,newFlushNum);


}



///
/// @param y
/// @param x center
/// @return known proposal function, which is normal distribution
double mc_computation::Q(const double &y, const double &x, const double &a, const double &b){

    double val=1/(std::pow(2.0*PI,0.5)*h)
               *std::exp(-1/(2*std::pow(h,2))*std::pow(y-x,2.0));

    return val;

}


///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return truncated Gaussian
double mc_computation::f(const double &y, const double &x, const double &a, const double &b){


    if(y<=a or y>=b){
        return 0;
    }else{

        double val=std::exp(-1.0/(2.0*std::pow(h,2))*std::pow(y-x,2));
        return val;
    }

}


///
/// @param x center
/// @param a left end
/// @param b right end
/// @return random number from truncated Gaussian
double mc_computation::reject_sampling_one_data(const double &x,const double &a, const double &b){

    std::random_device rd;  // Create a random device object
    std::ranlux24_base engine(rd());  // Seed the engine with the random device

    std::normal_distribution<> normal_dist(x,h);
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    double y=normal_dist(engine);
    double u=distUnif01(engine);

    while(u>=f(y,x,a,b)/(M* Q(y,x,a,b))){
        y=normal_dist(engine);
        u=distUnif01(engine);

    }

    return y;

}

///
/// @param x center
/// @param a left end
/// @param b right end
/// @return integral
double mc_computation::zVal(const double& x,const double &a, const double &b){

    auto integrandWithParam=[x,a,b, this](const double &y){
        return this->integrand(y,x,a,b);
    };
    double result = boost::math::quadrature::trapezoidal(integrandWithParam,a,b);

    return result;

}

///
/// @param y
/// @param x center
/// @param a left end
/// @param b right end
/// @return
double mc_computation::integrand(const double &y, const double& x,const double &a, const double &b){

    return f(y,x,a,b);


}

///
/// @param LCurr
/// @param d0VecCurr
/// @param d1VecCurr
/// @param LNext
/// @param d0VecNext
/// @param d1VecNext
/// @return
void mc_computation::proposal(const double &LCurr, const std::shared_ptr<double[]>& d0VecCurr ,const std::shared_ptr<double[]>&d1VecCurr,
              double & LNext, std::shared_ptr<double[]>& d0VecNext, std::shared_ptr<double[]>& d1VecNext){


    //proposal using truncated Gaussian
    double lm = potFuncPtr->getLm();
    for(int i=0;i<N;i++){
        d0VecNext[i]=reject_sampling_one_data(d0VecCurr[i],0,lm);
    }

    for(int i=0;i<N-1;i++){
        d1VecNext[i]=reject_sampling_one_data(d1VecCurr[i],0,lm);
    }

    LNext=reject_sampling_one_data(LCurr,0,lm);

}


///
/// @param LCurr
/// @param d0VecCurr
/// @param d1VecCurr
/// @param UCurr
/// @param LNext
/// @param d0VecNext
/// @param d1VecNext
/// @param UNext
/// @return
double mc_computation::acceptanceRatio(const double &LCurr,const std::shared_ptr<double[]>& d0VecCurr ,const std::shared_ptr<double[]>&d1VecCurr
        ,const double& UCurr,
                       const double &LNext, const std::shared_ptr<double[]>& d0VecNext,
                       const std::shared_ptr<double[]>&  d1VecNext,
                       double &UNext){

    double lm=potFuncPtr->getLm();

    UNext=(*potFuncPtr)(LNext,d0VecNext.get(),d1VecNext.get());

    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    double zLCurr= zVal(LCurr,0,lm);
    double zLNext= zVal(LNext,0,lm);

    double ratio_L=zLCurr/zLNext;
    R*=ratio_L;


    for(int i=0;i<N;i++){
    double zd0OneCurrVal= zVal(d0VecCurr[i],0,lm);
    double zd0OneNextVal= zVal(d0VecNext[i],0,lm);
    double ratio_d0OneVal=zd0OneCurrVal/zd0OneNextVal;
    R*=ratio_d0OneVal;

    }

    for(int i=0;i<N-1;i++){
        double zd1OneCurrVal= zVal(d1VecCurr[i],0,lm);
        double zd1OneNextVal= zVal(d1VecNext[i],0,lm);
        double ratio_d1OneVal=zd1OneCurrVal/zd1OneNextVal;
        R*=ratio_d1OneVal;
    }

    return std::min(1.0,R);



}
