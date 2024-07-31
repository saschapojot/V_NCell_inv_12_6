//
// Created by polya on 7/19/24.
//

#ifndef V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#define V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
#include "../potentialFunction/potentialFunctionPrototype.hpp"
#include <boost/filesystem.hpp>

#include <boost/math/quadrature/trapezoidal.hpp>
#include <chrono>
#include <cstdlib>
#include <cxxabi.h>
#include <fstream>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <typeinfo>
#include <vector>

namespace fs = boost::filesystem;
const auto PI=M_PI;


class mc_computation {
public:
    mc_computation(const std::string &cppInParamsFileName) {
        std::ifstream file(cppInParamsFileName);
        if (!file.is_open()) {
            std::cerr << "Failed to open the file." << std::endl;
            std::exit(20);
        }
        std::string line;

        int paramCounter = 0;

        while (std::getline(file, line)) {
            // Check if the line is empty
            if (line.empty()) {
                continue; // Skip empty lines
            }
            std::istringstream iss(line);

            //read T
            if (paramCounter == 0) {
                iss >> T;
                if (T <= 0) {
                    std::cerr << "T must be >0" << std::endl;
                    std::exit(1);
                }//end if
                std::cout << "T=" << T << std::endl;
                this->beta = 1 / T;
                double stepForT1 = 0.1;
                double h_threshhold=0.02;
                this->h=h_threshhold;
//                this->h = stepForT1 * T > h_threshhold ? h_threshhold : stepForT1 * T;//stepSize;
                std::cout << "h=" << h << std::endl;
                this->M = std::pow(2.0 * PI, 0.5) * h * 1.001;
                std::cout<<"M="<<M<<std::endl;
                paramCounter++;
                continue;
            }//end reading T
            //read coefficients
            if(paramCounter==1){
                iss>>coefsToPotFunc;
                paramCounter++;
                continue;

            }// end reading coefficients


            //read potential function name
            if(paramCounter==2){
                iss>>potFuncName;
                paramCounter++;
                continue;
            }//end reading potential function name

            //read initial values
            if(paramCounter==3){
                std::string temp;
                if (std::getline(iss, temp, ',')){
                    LInit=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    y0Init=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    z0Init=std::stod(temp);
                }
                if (std::getline(iss, temp, ',')){
                    y1Init=std::stod(temp);
                }
                paramCounter++;
                continue;




            }//end reading initial values


            //read loopToWrite
            if(paramCounter==4){
                iss>>loopToWrite;
                paramCounter++;
                continue;
            }//end reading loopToWrite

            //read newFlushNum
            if(paramCounter==5){
                iss>>newFlushNum;
                paramCounter++;
                continue;
            }//end reading newFlushNum

            //read loopLastFile
            if(paramCounter==6){
                //if loopLastFileStr is "-1", loopLastFile uses the overflowed value
                //and loopLastFile+1 will be 0
                iss>>loopLastFile;
                paramCounter++;
                continue;
            }//end reading loopLastFile

            //read TDirRoot
            if (paramCounter==7){
                iss>>TDirRoot;
                paramCounter++;
                continue;
            }//end reading TDirRoot

            //read U_dist_dataDir
            if(paramCounter==8){
                iss>>U_dist_dataDir;
                paramCounter++;
                continue;
            }//end reading U_dist_dataDir



        }//end while
        this->potFuncPtr = createPotentialFunction(potFuncName, coefsToPotFunc);
        potFuncPtr->init();
        this->varNum = 5;
        try {
            this->U_dist_ptr= std::shared_ptr<double[]>(new double[loopToWrite * varNum],
                                                        std::default_delete<double[]>());
        }
        catch (const std::bad_alloc &e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
            std::exit(2);
        } catch (const std::exception &e) {
            std::cerr << "Exception: " << e.what() << std::endl;
            std::exit(2);
        }
        std::cout<<"LInit="<<LInit<<", y0Init="<<y0Init
                 <<", z0Init="<<z0Init<<", y1Init="<<y1Init<<std::endl;

        std::cout<<"loopToWrite="<<loopToWrite<<std::endl;
        std::cout<<"newFlushNum="<<newFlushNum<<std::endl;
        std::cout<<"loopLastFile+1="<<loopLastFile+1<<std::endl;
        std::cout<<"TDirRoot="<<TDirRoot<<std::endl;
        std::cout<<"U_dist_dataDir="<<U_dist_dataDir<<std::endl;

    }//end constructor




public:
//    ///
//    /// @param x
//    /// @param leftEnd
//    /// @param rightEnd
//    /// @param eps
//    /// @return return a value within distance eps from x, on the open interval (leftEnd, rightEnd)
//   double generate_uni_open_interval(const double &x, const double &leftEnd, const double &rightEnd, const double &eps);
//


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
    bool proposal(const double &LCurr, const double& y0Curr,const double& z0Curr, const double& y1Curr,
                  double & LNext, double & y0Next, double & z0Next, double & y1Next, double &LReset);


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
    double acceptanceRatio(const double &LCurr,const double &y0Curr,
                           const double &z0Curr, const double& y1Curr,const double& UCurr,
                           const double &LNext, const double& y0Next,
                           const double & z0Next, const double & y1Next,
                           double &UNext,const double &LReset);
//    ///
//    /// @param x proposed value
//    /// @param y current value
//    /// @param a left end of interval
//    /// @param b right end of interval
//    /// @param epsilon half length
//    /// @return proposal probability S(x|y)
//    double S(const double &x, const double &y,const double &a, const double &b, const double &epsilon);

    ///
/// @param y
/// @param x center
/// @param a left end
///@param b right end
/// @return known proposal function, which is normal distribution
    double Q(const double &y, const double &x, const double &a, const double &b);

    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return truncated Gaussian
    double f(const double &y, const double &x, const double &a, const double &b);

    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return random number from truncated Gaussian
    double reject_sampling_one_data(const double &x,const double &a, const double &b);

    ///
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return integral
    double zVal(const double& x,const double &a, const double &b);

    ///
    /// @param y
    /// @param x center
    /// @param a left end
    /// @param b right end
    /// @return
    double integrand(const double &y, const double& x,const double &a, const double &b);

    void execute_mc(const double& L,const double &y0, const double &z0, const double& y1, const size_t & loopInit, const size_t & flushNum);


    static void saveArrayToCSV(const std::shared_ptr<double[]>& array, const  size_t& arraySize, const std::string& filename, const size_t& numbersPerRow) ;

    void init_and_run();

public:
    double T;// temperature
    double beta;
    double h;// step size
    size_t loopToWrite;
    size_t newFlushNum;
    size_t loopLastFile;
    std::shared_ptr<potentialFunction> potFuncPtr;
    std::string TDirRoot;
    std::string U_dist_dataDir;
    std::shared_ptr<double[]> U_dist_ptr;
    size_t varNum;
    double LInit;
    double y0Init;
    double z0Init;
    double y1Init;
    std::string coefsToPotFunc;
    std::string potFuncName;
    double M;

};


#endif //V_INV_12_6_MC_READ_LOAD_COMPUTE_HPP
