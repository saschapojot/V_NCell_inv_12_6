//
// Created by polya on 7/19/24.
//

#include "mc_read_load_compute.hpp"





void mc_computation::execute_mc(const double& L,const std::shared_ptr<double[]>& d1Vec, const std::shared_ptr<double[]>& d2Vec, const size_t & loopInit, const size_t & flushNum){

    double LCurr = L;
    std::shared_ptr<double[]> d1VecCurr=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

    std::shared_ptr<double[]> d2VecCurr=std::shared_ptr<double[]>(new double[N-1], std::default_delete<double[]>());

    std::memcpy(d1VecCurr.get(),d1Vec.get(),N*sizeof (double ));
    std::memcpy(d2VecCurr.get(),d2Vec.get(),(N-1)*sizeof (double ));

//    std::cout<<"d1VecCurr: ";
//    print_shared_ptr(d1VecCurr,N);
//    std::cout<<"d2VecCurr: ";
//    print_shared_ptr(d2VecCurr,N-1);

    //initialize next values
    double LNext;
    std::shared_ptr<double[]> d1VecNext=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());
    std::shared_ptr<double[]> d2VecNext=std::shared_ptr<double[]>(new double[N-1], std::default_delete<double[]>());

    double UCurr;
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

            this->proposal_unit(LCurr,d1VecCurr,d2VecCurr,LNext,d1VecNext,d2VecNext);
            double UNext;
            UCurr=(*potFuncPtr)(LCurr,d1VecCurr.get(),d2VecCurr.get());
            double r= acceptanceRatio_uni(LCurr,d1VecCurr,d2VecCurr,UCurr,
                                      LNext,d1VecNext,d2VecNext,UNext);
            double u = distUnif01(e2);
            if (u <= r) {
                LCurr = LNext;
                std::memcpy(d1VecCurr.get(),d1VecNext.get(),N*sizeof (double ));
                std::memcpy(d2VecCurr.get(),d2VecNext.get(),(N-1)*sizeof (double ));


                UCurr = UNext;

            }//end of accept-reject
            U_dist_ptr[varNum*j+0]=UCurr;
            U_dist_ptr[varNum*j+1]=LCurr;
            for(int n=2;n<=1+N;n++){
                U_dist_ptr[varNum*j+n]=d1VecCurr[n-2];
            }
            for(int n=N+2;n<=2*N;n++){
                U_dist_ptr[varNum*j+n]=d2VecCurr[n-N-2];
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
        outFile<<",d1"+std::to_string(i);
    }
    for(int i=0;i<N-1;i++){
        outFile<<",d2"+std::to_string(i);
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
    this->execute_mc(LInit,d1VecInit,d2VecInit,loopLastFile+1,newFlushNum);


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
/// @param d1VecCurr
/// @param d2VecCurr
/// @param LNext
/// @param d1VecNext
/// @param d2VecNext
/// @return
void mc_computation::proposal(const double &LCurr, const std::shared_ptr<double[]>& d1VecCurr ,const std::shared_ptr<double[]>&d2VecCurr,
              double & LNext, std::shared_ptr<double[]>& d1VecNext, std::shared_ptr<double[]>& d2VecNext){


    //proposal using truncated Gaussian
    double lm = potFuncPtr->getLm();
    for(int i=0;i<N;i++){
        d1VecNext[i]=reject_sampling_one_data(d1VecCurr[i],0,lm);
    }

    for(int i=0;i<N-1;i++){
        d2VecNext[i]=reject_sampling_one_data(d2VecCurr[i],0,lm);
    }

    LNext=reject_sampling_one_data(LCurr,0,lm);

}


///
/// @param LCurr
/// @param d1VecCurr
/// @param d2VecCurr
/// @param UCurr
/// @param LNext
/// @param d1VecNext
/// @param d2VecNext
/// @param UNext
/// @return
double mc_computation::acceptanceRatio(const double &LCurr,const std::shared_ptr<double[]>& d1VecCurr ,const std::shared_ptr<double[]>&d2VecCurr
        ,const double& UCurr,
                       const double &LNext, const std::shared_ptr<double[]>& d1VecNext,
                       const std::shared_ptr<double[]>&  d2VecNext,
                       double &UNext){

    double lm=potFuncPtr->getLm();

    UNext=(*potFuncPtr)(LNext,d1VecNext.get(),d2VecNext.get());

    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    double zLCurr= zVal(LCurr,0,lm);
    double zLNext= zVal(LNext,0,lm);

    double ratio_L=zLCurr/zLNext;
    R*=ratio_L;


    for(int i=0;i<N;i++){
    double zd1OneCurrVal= zVal(d1VecCurr[i],0,lm);
    double zd1OneNextVal= zVal(d1VecNext[i],0,lm);
    double ratio_d1OneVal=zd1OneCurrVal/zd1OneNextVal;
    R*=ratio_d1OneVal;

    }

    for(int i=0;i<N-1;i++){
        double zd2OneCurrVal= zVal(d2VecCurr[i],0,lm);
        double zd2OneNextVal= zVal(d2VecNext[i],0,lm);
        double ratio_d2OneVal=zd2OneCurrVal/zd2OneNextVal;
        R*=ratio_d2OneVal;
    }

    return std::min(1.0,R);



}


///
/// @param x proposed value
/// @param y current value
/// @param a left end of interval
/// @param b right end of interval
/// @param epsilon half length
/// @return proposal probability S(x|y)
double mc_computation::S_uni(const double &x, const double &y,const double &a, const double &b, const double &epsilon){

    if (a<y and y<a+epsilon){
        return 1.0/(y-a+epsilon);
    } else if( a+epsilon<=y and y<b+epsilon){
        return 1.0/(2.0*epsilon);
    }else if(b-epsilon<=y and y<b){
        return 1/(b-y+epsilon);
    } else{

        std::cerr<<"value out of range."<<std::endl;
        std::exit(10);


    }


}


///
/// @param LCurr
/// @param d1VecCurr
/// @param d2VecCurr
/// @param LNext
/// @param d1VecNext
/// @param d2VecNext
void mc_computation::proposal_unit(const double &LCurr, const std::shared_ptr<double[]>& d1VecCurr ,const std::shared_ptr<double[]>&d2VecCurr,
                   double & LNext, std::shared_ptr<double[]>& d1VecNext, std::shared_ptr<double[]>& d2VecNext){
    //proposal using uniform distribution
    double lm=potFuncPtr->getLm();
    for(int i=0;i<N;i++){
        d1VecNext[i]= generate_uni_open_interval(d1VecCurr[i],0,lm,h);
    }

    for(int i=0;i<N-1;i++){
        d2VecNext[i]=generate_uni_open_interval(d2VecCurr[i],0,lm,h);
    }

    LNext=generate_uni_open_interval(LCurr,0,lm,h);
}


///
/// @param x
/// @param leftEnd
/// @param rightEnd
/// @param eps
/// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps){


    double xMinusEps=x-eps;
    double xPlusEps=x+eps;

    double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
    double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;
//    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
//std::cout<<"x="<<x<<std::endl;
//std::cout<<"unif_left_end="<<unif_left_end<<std::endl;
//std::cout<<"unif_right_end="<<unif_right_end<<std::endl;
    std::random_device rd;
    std::ranlux24_base e2(rd());
// in std::uniform_real_distribution<> distUnif(a,b), the random numbers are from interval [a, b)
//we need random numbers from interval (a,b)
    double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());
//    std::cout<<"unif_left_end_double_on_the_right="<<unif_left_end_double_on_the_right<<std::endl;



    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)

    double xNext=distUnif(e2);
    return xNext;



}

///
/// @param LCurr
/// @param d1VecCurr
/// @param d2VecCurr
/// @param UCurr
/// @param LNext
/// @param d1VecNext
/// @param d2VecNext
/// @param UNext
/// @return
double mc_computation::acceptanceRatio_uni(const double &LCurr,const std::shared_ptr<double[]>& d1VecCurr ,const std::shared_ptr<double[]>&d2VecCurr
        ,const double& UCurr,
                           const double &LNext, const std::shared_ptr<double[]>& d1VecNext,
                           const std::shared_ptr<double[]>&  d2VecNext,
                           double &UNext) {


    double lm = potFuncPtr->getLm();
    UNext = (*potFuncPtr)(LNext, d1VecNext.get(), d2VecNext.get());
    double numerator = -this->beta * UNext;
    double denominator = -this->beta * UCurr;
    double R = std::exp(numerator - denominator);

    double S_LCurrNext = S_uni(LCurr, LNext, 0, lm, h);
    double S_LNextCurr = S_uni(LNext, LCurr, 0, lm, h);

    double ratio_L = S_LCurrNext / S_LNextCurr;
    if (std::fetestexcept(FE_DIVBYZERO)) {
        std::cout << "Division by zero exception caught." << std::endl;
        std::exit(15);
    }

    if (std::isnan(ratio_L)) {
        std::cout << "The result is NaN." << std::endl;
        std::exit(15);
    }
    R *= ratio_L;
    //end L ratio

    for(int i=0;i<N;i++){
        double S_d1iCurrNext= S_uni(d1VecCurr[i],d1VecNext[i],0,lm,h);
        double S_d1iNextCurr= S_uni(d1VecNext[i],d1VecCurr[i],0,lm,h);
        double ratio_d1i=S_d1iCurrNext/S_d1iNextCurr;
        if (std::fetestexcept(FE_DIVBYZERO)) {
            std::cout << "Division by zero exception caught." << std::endl;
            std::exit(15);
        }

        if (std::isnan(ratio_d1i)) {
            std::cout << "The result is NaN." << std::endl;
            std::exit(15);
        }
        R*=ratio_d1i;



    }//end d1 ratio

    for(int i=0;i<N-1;i++){
        double S_d2iCurrNext= S_uni(d2VecCurr[i],d2VecNext[i],0,lm,h);
        double S_d2iNextCurr= S_uni(d2VecNext[i],d2VecCurr[i],0,lm,h);
        double ratio_d2i=S_d2iCurrNext/S_d2iNextCurr;
        if (std::fetestexcept(FE_DIVBYZERO)) {
            std::cout << "Division by zero exception caught." << std::endl;
            std::exit(15);
        }

        if (std::isnan(ratio_d2i)) {
            std::cout << "The result is NaN." << std::endl;
            std::exit(15);
        }
        R*=ratio_d2i;

    }//end ratio d2
    return std::min(1.0,R);



}