#ifndef AUXIL_HPP_
#define AUXIL_HPP_

#include <fstream>
#include <iostream>
#include <time.h>
#include <random>
#include <array> 


#include <baobzi_template.hpp>

struct timespec get_wtime() {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts;
}

double get_wtime_diff(const struct timespec *ts, const struct timespec *tf) {
    return (tf->tv_sec - ts->tv_sec) + (tf->tv_nsec - ts->tv_nsec) * 1E-9;
}

double mean_error(std::vector<double> const& xarray) {
    auto const count = static_cast<double>(xarray.size());
    return std::reduce(xarray.begin(), xarray.end()) / count; 
}

double total_sum(std::vector<double> const& xarray) {
    auto const count = static_cast<double>(xarray.size());
    return std::reduce(xarray.begin(), xarray.end()); 
}

void speak(const char name[], const double var) {
    std::cout << name << ": " << var << std::endl;
    return; 
}

void speakvec(std::vector<double> thevec) {
    size_t n = thevec.size();
    for (size_t i = 0; i < n; i++) std::cout << thevec[i] << ","; 
    std::cout << std::endl; 
}

double find_order(const double input) {
    return powf(10.0f, floorf(log10f(input)));
}

std::vector<double> createErr(const int a, const int b) {
    std::vector<double> errVec; 
    for (int i = b; i >= a; i--) {
        errVec.push_back(10 ^ (-i)); 
    }
    return errVec; 
}
int binarySearch(const std::vector<double>& numbers, int target, int start, int end) {
    if (start > end) {
        return -1;
    }

    int mid = start + (end - start) / 2;

    if (numbers[mid] <= target && target <= numbers[mid + 1]) {
        return mid;
    } 
    else if (numbers[mid] < target) {
        return binarySearch(numbers, target, mid + 1, end);
    } 
    else {
        return binarySearch(numbers, target, start, mid - 1);
    }
}

int findIntervalIndex(const std::vector<double>& numbers, int target) {
    if (target < numbers[0] || target > numbers[numbers.size() - 1]) {
        return -1;
    }

    return binarySearch(numbers, target, 0, numbers.size() - 1);
}
















#endif
