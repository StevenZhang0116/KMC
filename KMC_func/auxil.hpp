/**
 * @file auxil.hpp
 * @author Zihan Zhang
 * @brief Auxiliary functions that are used in construct Baobzi/Baobzi family object; 
 *        the purposes for most functions are self-evident by naming
 *
 * @copyright Copyright (c) 2023
 *
 */

#ifndef AUXIL_HPP_
#define AUXIL_HPP_

#include <fstream>
#include <iostream>
#include <time.h>
#include <random>
#include <array> 
#include <numeric>
#include <string>
#include <filesystem>


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

template <typename T>
void speak(const char name[], const T& var) {
    std::cout << name << ": " << var << std::endl;
}

template <typename T>
void printElement(const T& element) {
    std::cout << element << " ";
}

template <typename T>
void printVector(const std::vector<T>& vec) {
    for (const auto& element : vec) {
        printElement(element);
    }
    std::cout << std::endl;
}

template <typename T>
std::vector<T> getUniqueElements(const std::vector<T>& vec) {
    std::set<T> uniqueSet(vec.begin(), vec.end());
    std::vector<T> uniqueVec(uniqueSet.begin(), uniqueSet.end());
    return uniqueVec;
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
    if (start > end) return -1;

    int mid = start + (end - start) / 2;

    if (numbers[mid] <= target && target <= numbers[mid + 1]) return mid;
    else if (numbers[mid] < target) return binarySearch(numbers, target, mid + 1, end);
    else return binarySearch(numbers, target, start, mid - 1);
}

int findIntervalIndex(const std::vector<double>& numbers, int target) {
    if (target < numbers[0] || target > numbers[numbers.size() - 1]) return -1;

    return binarySearch(numbers, target, 0, numbers.size() - 1);
}

// find the minimum [0] and maximum [1] of a 2D matrix 
std::vector<double> findMinMaxVec(const std::vector<std::vector<double>> matrix) {
    double largestValue = 0;
    double smallestValue = 0; 

    for (const auto& row : matrix) {
        for (const auto& element : row) {
            if (element > largestValue) {
                largestValue = element;
            }
        }
    }

    for (const auto& row : matrix) {
        for (const auto& element : row) {
            if (element > largestValue) {
                smallestValue = element;
            }
        }
    }

    std::vector<double> res = {smallestValue, largestValue}; 
    return res; 
}

int round10(int n) {
    int a = (n / 10) * 10;
    int b = a + 10;  
    return (n - a > b - n)? b : a;
}

std::vector<double> cumulativeSum(const std::vector<double>& input) {
    std::vector<double> cumulative;
    double sum = 0;

    for (const auto& element : input) {
        sum += element;
        cumulative.push_back(sum);
    }
    return cumulative;
}

std::vector<double> readDataFromFile(const std::string& filename) {
    std::vector<double> data;
    std::ifstream inputFile(filename);

    if (inputFile.is_open()) {
        std::string line;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);
            std::string token;
            while (std::getline(iss, token, ',')) {
                // speak("token", token); 
                std::locale originalLocale(std::locale::classic());
                std::locale::global(std::locale::classic());
                double value = std::stod(token);
                data.push_back(value);
                std::locale::global(originalLocale);
            }
        }

        inputFile.close();
    } 
    else std::cout << "Failed to open the file." << std::endl;
    return data;
}

std::uintmax_t calculateFolderSize(const std::string& folderPath) {
    std::uintmax_t totalSize = 0;

    for (const auto& entry : std::filesystem::recursive_directory_iterator(folderPath)) {
        if (std::filesystem::is_regular_file(entry)) {
            totalSize += std::filesystem::file_size(entry);
        }
    }
    
    return totalSize;
}

double round_up(double value, int decimal_places) {
    const double multiplier = std::pow(10.0, decimal_places);
    return std::ceil(value * multiplier) / multiplier;
}


void removeAllFilesInFolder(const std::string& folderPath) {
    try {
        for (const auto& entry : std::filesystem::directory_iterator(folderPath)) {
            if (std::filesystem::is_regular_file(entry)) {
                std::filesystem::remove(entry.path());
            }
        }
    } catch (const std::filesystem::filesystem_error& ex) {
        std::cerr << "Error while removing files: " << ex.what() << std::endl;
    }
}




#endif
