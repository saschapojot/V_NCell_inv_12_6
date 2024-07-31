//
// Created by polya on 7/19/24.
//

#include "mc_read_load_compute.hpp"

/////
///// @param x
///// @param leftEnd
///// @param rightEnd
///// @param eps
///// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
//double mc_computation::generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps){
//
//
//double xMinusEps=x-eps;
//double xPlusEps=x+eps;
//
//double unif_left_end=xMinusEps<leftEnd?leftEnd:xMinusEps;
//double unif_right_end=xPlusEps>rightEnd?rightEnd:xPlusEps;
////    std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
////std::cout<<"x="<<x<<std::endl;
////std::cout<<"unif_left_end="<<unif_left_end<<std::endl;
////std::cout<<"unif_right_end="<<unif_right_end<<std::endl;
//    std::random_device rd;
//    std::ranlux24_base e2(rd());
//// in std::uniform_real_distribution<> distUnif(a,b), the random numbers are from interval [a, b)
////we need random numbers from interval (a,b)
//double unif_left_end_double_on_the_right=std::nextafter(unif_left_end, std::numeric_limits<double>::infinity());
////    std::cout<<"unif_left_end_double_on_the_right="<<unif_left_end_double_on_the_right<<std::endl;
//
//
//
//    std::uniform_real_distribution<> distUnif(unif_left_end_double_on_the_right,unif_right_end); //[unif_left_end_double_on_the_right, unif_right_end)
//
//    double xNext=distUnif(e2);
//    return xNext;
//
//
//
//}

/////
///// @param LCurr current value of L
///// @param y0Curr current value of y0
///// @param z0Curr current value of z0
///// @param y1Curr current value of y1
///// @param LNext  next value of L
///// @param y0Next next value of y0
///// @param z0Next next value of z0
///// @param y1Next next value of y1
//bool mc_computation::proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
//              double & LNext, double & y0Next, double & z0Next, double & y1Next, double &LReset){
//
//    double eps=potFuncPtr->get_eps();
//    double lm=potFuncPtr->getLm();
//    bool reset= false;
//
//    y0Next= generate_uni_open_interval(y0Curr,0,LCurr,eps);
//    z0Next= generate_uni_open_interval(z0Curr,0,LCurr,eps);
//    y1Next= generate_uni_open_interval(y1Curr,0,LCurr,eps);
//    //to prevent the case that LCurr<=y0Next+z0Next+y1Next,
//    //or the case that LCurr>=lm,
//    // we set LCurr in (y0Next+z0Next+y1Next, lm)
//    double LCurrReset=LCurr;
//    if(LCurr<=y0Next+z0Next+y1Next or LCurr>=lm) {
//        std::random_device rd;
//        std::ranlux24_base e2(rd());
//        double past_left = std::nextafter(y0Next + z0Next + y1Next, std::numeric_limits<double>::infinity());
//        std::uniform_real_distribution<> distUnif(past_left, lm);
//        LCurrReset = distUnif(e2);
//        reset= true;
//
//        LReset=LCurrReset;
//
//    }
//
//    LNext= generate_uni_open_interval(LCurrReset,y0Next+z0Next+y1Next,lm,eps);
//    LReset=LCurrReset;
//    return reset;
//
//
//
//
//}


//////
///// @param LCurr
///// @param y0Curr
///// @param z0Curr
///// @param y1Curr
///// @param UCurr
///// @param LNext
///// @param y0Next
///// @param z0Next
///// @param y1Next
///// @param UNext
///// @param LReset
///// @return
//double mc_computation::acceptanceRatio(const double &LCurr,const double &y0Curr,
//                       const double &z0Curr, const double& y1Curr,const double& UCurr,
//                       const double &LNext, const double& y0Next,
//                       const double & z0Next, const double & y1Next,
//                       double &UNext,const double &LReset){
//    double eps=potFuncPtr->get_eps();
//    double lm=potFuncPtr->getLm();
//
//    UNext=((*potFuncPtr)(LNext,y0Next,z0Next,y1Next));
//    double numerator = -this->beta*UNext;
//
//    double denominator=-this->beta*UCurr;
//
//    double R=std::exp(numerator - denominator);
//
//    double ratioL= S(LReset,LNext,y0Next+z0Next+y1Next,lm,eps)/S(LNext,LReset,y0Next+z0Next+y1Next,lm,eps);
//
//    double ratio_y0=S(y0Curr,y0Next,0,LCurr,eps)/S(y0Next,y0Curr,0,LCurr,eps);
//
//    double ratio_z0=S(z0Curr,z0Next,0,LCurr,eps)/S(z0Next,z0Curr,0,LCurr,eps);
//
//    double  ratio_y1=S(y1Curr,y1Next,0,LCurr,eps)/S(y1Next,y1Curr,0,LCurr,eps);
//
//    R*=ratioL*ratio_y0*ratio_z0*ratio_y1;
//
//
//    return std::min(1.0,R);
//
//
//
//
//}


/////
///// @param x proposed value
///// @param y current value
///// @param a left end of interval
///// @param b right end of interval
///// @param epsilon half length
///// @return proposal probability S(x|y)
//double mc_computation::S(const double &x, const double &y,const double &a, const double &b, const double &epsilon){
//
//    if (a<y and y<a+epsilon){
//        return 1.0/(y-a+epsilon);
//    } else if( a+epsilon<=y and y<b+epsilon){
//        return 1.0/(2.0*epsilon);
//    }else if(b-epsilon<=y and y<b){
//        return 1/(b-y+epsilon);
//    } else{
//
//        std::cerr<<"value out of range."<<std::endl;
//        std::exit(10);
//
//
//    }
//
//
//}
//


void mc_computation::execute_mc(const double& L,const double &y0, const double &z0, const double& y1, const size_t & loopInit, const size_t & flushNum){

    double LCurr = L;
    double y0Curr = y0;
    double z0Curr = z0;
    double y1Curr = y1;
    std::cout<<"Before mc: "<<"LCurr="<<LCurr<<", y0Curr="<<y0Curr<<", z0Curr="<<z0Curr<<", y1Curr="<<y1Curr<<std::endl;
    double UCurr;// = (*potFuncPtr)(LCurr, y0Curr, z0Curr, y1Curr);
    std::random_device rd;
    std::ranlux24_base e2(rd());
    std::uniform_real_distribution<> distUnif01(0, 1);//[0,1)
    size_t loopStart = loopInit;
    for (size_t fls = 0; fls < flushNum; fls++) {
        const auto tMCStart{std::chrono::steady_clock::now()};
        for (size_t j = 0; j < loopToWrite; j++) {
            //propose a move
            double LNext;
            double y0Next;
            double z0Next;
            double y1Next;
            double LReset;

            this->proposal(LCurr,y0Curr,z0Curr,y1Curr,LNext,y0Next,z0Next,y1Next,LReset);
            double UNext;
            UCurr=((*potFuncPtr))(LReset,y0Curr,z0Curr,y1Curr);
            double r= acceptanceRatio(LCurr,y0Curr,z0Curr,y1Curr,UCurr,LNext,y0Next,z0Next,y1Next,UNext,LReset);
            double u = distUnif01(e2);
            if (u <= r) {
                LCurr = LNext;
                y0Curr = y0Next;
                z0Curr = z0Next;
                y1Curr = y1Next;
                UCurr = UNext;

            }//end of accept-reject
            U_dist_ptr[varNum*j+0]=UCurr;
            U_dist_ptr[varNum*j+1]=LCurr;
            U_dist_ptr[varNum*j+2]=y0Curr;
            U_dist_ptr[varNum*j+3]=z0Curr;
            U_dist_ptr[varNum*j+4]=y1Curr;
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







void mc_computation::saveArrayToCSV(const std::shared_ptr<double[]>& array, const  size_t& arraySize, const std::string& filename, const size_t& numbersPerRow) {

    std::ofstream outFile(filename);

    if (!outFile.is_open()) {
        std::cerr << "Error opening file: " << filename << std::endl;
        return;
    }
    outFile << std::setprecision(std::numeric_limits<double>::digits10 + 1) << std::fixed;

    outFile<<"U,"<<"L,"<<"y0,"<<"z0,"<<"y1"<<"\n";
    for (size_t i = 0; i < arraySize; ++i) {
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
    this->execute_mc(LInit,y0Init,z0Init,y1Init,loopLastFile+1,newFlushNum);


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
/// @param LCurr current value of L
/// @param y0Curr current value of y0
/// @param z0Curr current value of z0
/// @param y1Curr current value of y1
/// @param LNext  next value of L
/// @param y0Next next value of y0
/// @param z0Next next value of z0
/// @param y1Next next value of y1
/// @param LReset current value resetted
bool mc_computation::proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
              double & LNext, double & y0Next, double & z0Next, double & y1Next, double &LReset) {

    //proposal using truncated Gaussian
    double lm = potFuncPtr->getLm();
    bool reset = false;
    y0Next = reject_sampling_one_data(y0Curr, 0, lm);
    z0Next = reject_sampling_one_data(z0Curr, 0, lm);
    y1Next = reject_sampling_one_data(y1Curr, 0, lm);
    //to prevent the case that LCurr<=y0Next+z0Next+y1Next,
    //or the case that LCurr>=lm,
    // we set LCurr in (y0Next+z0Next+y1Next, lm)
    double LCurrReset = LCurr;
//    if (LCurr <= y0Next + z0Next + y1Next or LCurr >= lm) {
//        std::random_device rd;
//        std::ranlux24_base e2(rd());
//        double past_left = std::nextafter(y0Next + z0Next + y1Next, std::numeric_limits<double>::infinity());
//        std::uniform_real_distribution<> distUnif(past_left, lm);
//        LCurrReset = distUnif(e2);
//        reset = true;
//
//        LReset = LCurrReset;
//
//    }
    LReset = LCurrReset;
    LNext = reject_sampling_one_data(LCurr, 0, lm);

    return reset;

}

///
/// @param LCurr
/// @param y0Curr
/// @param z0Curr
/// @param y1Curr
/// @param UCurr
/// @param LNext
/// @param y0Next
/// @param z0Next
/// @param y1Next
/// @param UNext
/// @param LReset
/// @return
double mc_computation::acceptanceRatio(const double &LCurr,const double &y0Curr,
                       const double &z0Curr, const double& y1Curr,const double& UCurr,
                       const double &LNext, const double& y0Next,
                       const double & z0Next, const double & y1Next,
                       double &UNext,const double &LReset){



    double lm=potFuncPtr->getLm();
    UNext=((*potFuncPtr)(LNext,y0Next,z0Next,y1Next));
    double numerator = -this->beta*UNext;
    double denominator=-this->beta*UCurr;
    double R=std::exp(numerator - denominator);

    double zLCurr= zVal(LCurr,0,lm);
    double zLNext= zVal(LNext,0,lm);

    double ratio_L=zLCurr/zLNext;


    double zy0Curr= zVal(y0Curr,0,lm);
    double zy0Next= zVal(y0Next,0,lm);

    double ratio_y0=zy0Curr/zy0Next;

    double zz0Curr= zVal(z0Curr,0,lm);
    double zz0Next= zVal(z0Next,0,lm);

    double ratio_z0=zz0Curr/zz0Next;

    double zy1Curr= zVal(y1Curr,0,lm);
    double zy1Next= zVal(y0Next,0,lm);

    double ratio_y1=zy1Curr/zy1Next;

    R*=ratio_L*ratio_y0*ratio_z0*ratio_y1;

    return std::min(1.0,R);

}