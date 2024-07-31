#include <immintrin.h>
#include <iostream>
#include <memory>
#include <chrono>
#include <cmath>


// Custom approximate vectorized log function using AVX2
__m256d custom_mm256_log_pd(__m256d x) {
    alignas(32) double temp[4];
    _mm256_storeu_pd(temp, x);
    for (int i = 0; i < 4; ++i) {
        temp[i] = std::log(temp[i]);  // Using std::log for demonstration purposes
    }
    return _mm256_loadu_pd(temp);
}

// Custom approximate vectorized exp function using AVX2
__m256d custom_mm256_exp_pd(__m256d x) {
    alignas(32) double temp[4];
    _mm256_storeu_pd(temp, x);
    for (int i = 0; i < 4; ++i) {
        temp[i] = std::exp(temp[i]);  // Using std::exp for demonstration purposes
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
// Plain for loop power function
void plain_pow(const double* base, double exp, double* result, int n) {
    for (int i = 0; i < n; ++i) {
        result[i] = std::pow(base[i], exp);
    }
}

int main(int argc, char *argv[]) {
    int n = 100;  // Example size
    double exponent = 3.0;
    std::shared_ptr<double[]> base(new double[n]);
    std::shared_ptr<double[]> result_plain(new double[n]);
    std::shared_ptr<double[]> result_vectorized(new double[n]);

    // Initialize the base array with some values
    for (int i = 0; i < n; ++i) {
        base[i] = static_cast<double>(i + 1);
    }

    // Measure time for plain power function
    auto start = std::chrono::high_resolution_clock::now();
    plain_pow(base.get(), exponent, result_plain.get(), n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> plain_time = end - start;

    // Measure time for vectorized power function
    start = std::chrono::high_resolution_clock::now();
    vectorized_pow(base.get(), exponent, result_vectorized.get(), n);
    end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> vectorized_time = end - start;

    // Output results and timings
    std::cout << "Plain power function time: " << plain_time.count() << " seconds\n";
    std::cout << "Vectorized power function time: " << vectorized_time.count() << " seconds\n";

    // Compute and output the speedup
    double speedup = plain_time.count() / vectorized_time.count();
    std::cout << "Speedup: " << speedup << "x\n";

    // Verify the result for the first few elements
    bool correct = true;
    for (int i = 0; i < 10; ++i) {
        if (std::abs(result_plain[i] - result_vectorized[i]) > 1e-9) {
            correct = false;
            break;
        }
    }
    if (correct) {
        std::cout << "Results are correct for the first 10 elements.\n";
    } else {
        std::cout << "Results differ for the first 10 elements.\n";
    }

    return 0;
}