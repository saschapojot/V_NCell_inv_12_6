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

    // Vectorized multiplication summation function
    double vectorized_multiply_sum(const double* array, double multiplier, int n) {
        __m256d vec_multiplier = _mm256_set1_pd(multiplier); // Set all elements of the vector to the multiplier
        __m256d vec_sum = _mm256_setzero_pd(); // Initialize a 256-bit vector with zeros
        int i;

        // Process the array in chunks of 4 elements
        for (i = 0; i <= n - 4; i += 4) {
            __m256d vec_array = _mm256_loadu_pd(&array[i]); // Load 4 elements from the array
            __m256d vec_result = _mm256_mul_pd(vec_array, vec_multiplier); // Multiply each element by the multiplier
            vec_sum = _mm256_add_pd(vec_sum, vec_result); // Accumulate the sum
        }

        // Handle remaining elements
        double sum = 0.0;
        for (; i < n; ++i) {
            sum += array[i] * multiplier;
        }

        // Horizontal addition of the elements in vec_sum
        alignas(32) double temp[4];
        _mm256_storeu_pd(temp, vec_sum);
        for (int j = 0; j < 4; ++j) {
            sum += temp[j];
        }

        return sum;
    }

    // Custom approximate vectorized log function using AVX2
    __m256d custom_mm256_log_pd(__m256d x) {
        alignas(32) double temp[4];
        _mm256_storeu_pd(temp, x);
        for (int i = 0; i < 4; ++i) {
            temp[i] = std::log(temp[i]);
        }
        return _mm256_loadu_pd(temp);
    }

    // Custom approximate vectorized exp function using AVX2
    __m256d custom_mm256_exp_pd(__m256d x) {
        alignas(32) double temp[4];
        _mm256_storeu_pd(temp, x);
        for (int i = 0; i < 4; ++i) {
            temp[i] = std::exp(temp[i]);
        }
        return _mm256_loadu_pd(temp);
    }

    // Custom vectorized power function using AVX2
    __m256d custom_mm256_pow_pd(__m256d base, __m256d exp) {
        __m256d log_base = custom_mm256_log_pd(base); // Compute log(base)
        __m256d log_base_times_exp = _mm256_mul_pd(log_base, exp); // Compute log(base) * exp
        return custom_mm256_exp_pd(log_base_times_exp); // Compute exp(log(base) * exp)
    }

// Vectorized power function using a single double value as exponent
    void vectorized_pow(const double* base, double exp, double* result, int n) {
        __m256d vec_exp = _mm256_set1_pd(exp); // Set all elements of the vector to the exponent value
        int i;

        // Process the array in chunks of 4 elements
        for (i = 0; i <= n - 4; i += 4) {
            __m256d base_vals = _mm256_loadu_pd(&base[i]);
            __m256d result_vals = custom_mm256_pow_pd(base_vals, vec_exp);
            _mm256_storeu_pd(&result[i], result_vals);
        }

        // Process remaining elements
        for (; i < n; ++i) {
            result[i] = std::pow(base[i], exp);
        }
    }

    void init() override{
        this->json2Coefs(coefsInStr);
        this->r1=std::pow(2.0*a1/b1,1.0/6.0);
        this->r2=std::pow(2.0*a2/b2,1.0/6.0);
        this->lm=(static_cast<double >(N)*(r1+r2))*1.5;
        this->eps=((r1+r2)/2.0)/8;
        pow_result_tmp=std::shared_ptr<double[]>(new double[N], std::default_delete<double[]>());

        std::cout << "a1=" << a1 << ", b1=" << b1 << ", a2=" << a2 << ", b2=" << b2 << std::endl;
        std::cout<<"r1="<<r1<<", r2="<<r2<<std::endl;
        std::cout<<"lm="<<lm<<std::endl;
//        std::cout<<"eps="<<eps<<std::endl;


    }

    ///
    /// @param d0Vec array of distance d0, size=N
    /// @param d1Vec array of distance d1, size=N
    /// @param L total length
    /// @param N unit cell number
    /// @return
    double operator()(const double *d0Vec, double *d1Vec, const double&L) override {
//////////////////////////vectorization
//    //the last element of d1Vec is compute from d0Vec, d1Vec's first N-1 elements, and L
//    double d0SumNeg=vectorized_multiply_sum(d0Vec,-1,N);
////    std::cout<<"d0SumNeg="<<d0SumNeg<<std::endl;
//    double d1SumNeg= vectorized_multiply_sum(d1Vec,-1,N-1);
////    std::cout<<"d1SumNeg="<<d1SumNeg<<std::endl;
//    d1Vec[N-1]=d0SumNeg+d1SumNeg+L;
////    std::cout<<"d1Vec[N-1]="<<d1Vec[N-1]<<std::endl;
//
//
//    double val=0;
//    //d0Vec, -12-th power
//        vectorized_pow(d0Vec,-12.0,pow_result_tmp.get(),N);
//        val+= vectorized_multiply_sum(pow_result_tmp.get(),a1,N);
//        //d0Vec, -6-th power
//        vectorized_pow(d0Vec,-6.0,pow_result_tmp.get(),N);
//        val+=vectorized_multiply_sum(pow_result_tmp.get(),-b1,N);
//
//        //d1Vec, -12-th power
//        vectorized_pow(d1Vec,-12.0,pow_result_tmp.get(),N);
//        val+= vectorized_multiply_sum(pow_result_tmp.get(),a2,N);
//
//        //d1Vec, -6th power
//        vectorized_pow(d1Vec,-6.0,pow_result_tmp.get(),N);
//        val+=vectorized_multiply_sum(pow_result_tmp.get(),-b2,N);
//////////////////////////end of vectorization
        double sum=0;
        for(int i=0;i<N;i++){
            sum+=-d0Vec[i];
        }
        for(int i=0;i<N-1;i++){
            sum+=-d1Vec[i];
        }
        sum+=L;
        d1Vec[N-1]=sum;

        double val=0;
        for(int i=0;i<N;i++){
            val+=a1*std::pow(d0Vec[i],-12);
        }
        for(int i=0;i<N;i++){

            val+=-b1*std::pow(d0Vec[i],-6);
        }

        for(int i=0;i<N;i++){
            val+=a2*std::pow(d1Vec[i],-12);
        }
        for(int i=0;i<N;i++){
            val+=-b2*std::pow(d1Vec[i],-6);
        }
        return val;


    }//end of () operator

    double plain_for(const double *d0Vec, double *d1Vec, const double&L)override{
        double sum=0;
        for(int i=0;i<N;i++){
            sum+=-d0Vec[i];
        }
        for(int i=0;i<N-1;i++){
            sum+=-d1Vec[i];
        }
        sum+=L;
        d1Vec[N-1]=sum;

        double val=0;
        for(int i=0;i<N;i++){
            val+=a1*std::pow(d0Vec[i],-12);
        }
        for(int i=0;i<N;i++){

            val+=-b1*std::pow(d0Vec[i],-6);
        }

        for(int i=0;i<N;i++){
            val+=a2*std::pow(d1Vec[i],-12);
        }
        for(int i=0;i<N;i++){
            val+=-b2*std::pow(d1Vec[i],-6);
        }
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