//
// Created by polya on 7/19/24.
//

#include "potentialFunctionPrototype.hpp"

class V_inv_12_6:public potentialFunction {

public:
    V_inv_12_6(const std::string &coefsStr):potentialFunction(){
        this->coefsInStr=coefsStr;
    }

public:
    void json2Coefs(const std::string &coefsStr)override{
        std::stringstream iss;
        iss<<coefsStr;
        std::string temp;
        //read a1
        if (std::getline(iss, temp, ',')){
            this->a1=std::stod(temp);
        }

        //read b1

        if (std::getline(iss, temp, ',')){
            this->b1=std::stod(temp);
        }

        //read a2
        if (std::getline(iss, temp, ',')){
            this->a2=std::stod(temp);
        }

        //read b2

        if (std::getline(iss, temp, ',')){
            this->b2=std::stod(temp);
        }

        //read N

        if (std::getline(iss, temp, ',')){
            this->N=std::stoi(temp);
        }
    }


    void init() override{
        this->json2Coefs(coefsInStr);
        this->r1=std::pow(2.0*a1/b1,1.0/6.0);
        this->r2=std::pow(2.0*a2/b2,1.0/6.0);
        this->lm=(static_cast<double >(N)*(r1+r2))*1.5;
        this->eps=((r1+r2)/2.0)/8;
//        pow_result_tmp=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

        std::cout << "a1=" << a1 << ", b1=" << b1 << ", a2=" << a2 << ", b2=" << b2 << std::endl;
        std::cout<<"r1="<<r1<<", r2="<<r2<<std::endl;
        std::cout<<"lm="<<lm<<std::endl;
//        std::cout<<"eps="<<eps<<std::endl;


    }

    ///
    /// @param d1Vec array of distance d1, size=N
    /// @param d2Vec array of distance d2, size=N
    /// @param L total length
    /// @param N unit cell number
    /// @return
    double operator()(const double&L,const double *d1Vec,const double *d2Vec) override {
//////////////////////////vectorization
//    //the last element of d2Vec is compute from d1Vec, d2Vec's first N-1 elements, and L
//    double d1SumNeg=vectorized_multiply_sum(d1Vec,-1,N);
////    std::cout<<"d1SumNeg="<<d1SumNeg<<std::endl;
//    double d2SumNeg= vectorized_multiply_sum(d2Vec,-1,N-1);
////    std::cout<<"d2SumNeg="<<d2SumNeg<<std::endl;
//    d2Vec[N-1]=d1SumNeg+d2SumNeg+L;
////    std::cout<<"d2Vec[N-1]="<<d2Vec[N-1]<<std::endl;
//
//
//    double val=0;
//    //d1Vec, -12-th power
//        vectorized_pow(d1Vec,-12.0,pow_result_tmp.get(),N);
//        val+= vectorized_multiply_sum(pow_result_tmp.get(),a1,N);
//        //d1Vec, -6-th power
//        vectorized_pow(d1Vec,-6.0,pow_result_tmp.get(),N);
//        val+=vectorized_multiply_sum(pow_result_tmp.get(),-b1,N);
//
//        //d2Vec, -12-th power
//        vectorized_pow(d2Vec,-12.0,pow_result_tmp.get(),N);
//        val+= vectorized_multiply_sum(pow_result_tmp.get(),a2,N);
//
//        //d2Vec, -6th power
//        vectorized_pow(d2Vec,-6.0,pow_result_tmp.get(),N);
//        val+=vectorized_multiply_sum(pow_result_tmp.get(),-b2,N);
//////////////////////////end of vectorization
        double sum=0;
        for(int i=0;i<N;i++){
            sum+=-d1Vec[i];
        }
        for(int i=0;i<N-1;i++){
            sum+=-d2Vec[i];
        }
        sum+=L;


        double val=0;
        for(int i=0;i<N;i++){
            val+=a1*std::pow(d1Vec[i],-12);
        }
        for(int i=0;i<N;i++){

            val+=-b1*std::pow(d1Vec[i],-6);
        }

        for(int i=0;i<N-1;i++){
            val+=a2*std::pow(d2Vec[i],-12);
        }
        for(int i=0;i<N-1;i++){
            val+=-b2*std::pow(d2Vec[i],-6);
        }
        val+=a2*std::pow(sum,-12)-b2*std::pow(sum,-6);
        return val;


    }//end of () operator

    double plain_for(const double&L,const double *d1Vec,const double *d2Vec)override{
        double sum=0;
        for(int i=0;i<N;i++){
            sum+=-d1Vec[i];
        }
        for(int i=0;i<N-1;i++){
            sum+=-d2Vec[i];
        }
        sum+=L;


        double val=0;
        for(int i=0;i<N;i++){
            val+=a1*std::pow(d1Vec[i],-12);
        }
        for(int i=0;i<N;i++){

            val+=-b1*std::pow(d1Vec[i],-6);
        }

        for(int i=0;i<N-1;i++){
            val+=a2*std::pow(d2Vec[i],-12);
        }
        for(int i=0;i<N-1;i++){
            val+=-b2*std::pow(d2Vec[i],-6);
        }
        val+=a2*std::pow(sum,-12)-b2*std::pow(sum,-6);
        return val;

    }

    double getLm() const override {
        return lm;
    }
    double get_eps() const override {
        return eps;
    }
public:
    double a1;
    double a2;
    double b1;
    double b2;
    std::string coefsInStr;
    double r1;//min position of V1
    double r2;//min position of V2
    double lm;//range of distances
    double eps;//half interval length of uniform distribution
    int N;
    std::shared_ptr<double[]>pow_result_tmp;
};

std::shared_ptr<potentialFunction>  createPotentialFunction(const std::string& funcName, const std::string &coefsJsonStr) {
    if (funcName == "V_inv_12_6") {

        return std::make_shared<V_inv_12_6>(coefsJsonStr);
    }

    else {
        throw std::invalid_argument("Unknown potential function type");
    }
}