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

double mean_error(const double xarray[], int gridcnt) {
    return std::accumulate(xarray, xarray + gridcnt, 0) / gridcnt;
}

void speak(const char name[], const double var) {
    std::cout << name << ", " << var << std::endl;
    return; 
}

void speakvec(const double varvec[], size_t len) {
    for (size_t i = 0; i < len; i++){
        std::cout << varvec[i] << ","; 
    }
    std::cout << std::endl; 
    std::cout << "==============" << std::endl; 
}

#endif
