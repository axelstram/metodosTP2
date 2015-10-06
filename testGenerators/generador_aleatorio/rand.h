#ifndef RAND_H
#define RAND_H

#include <iostream>
#include <functional>
#include <random>
#include <chrono>
#include <stdio.h>
#include <stdlib.h>


class UniformDist {
    public:

        UniformDist(double min, double max) : distribution(min, max), r(std::bind(distribution, std::default_random_engine(seed))) {}
        double Generate() { return r(); }

    private:

        std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::uniform_real_distribution<double> distribution;
        std::function<double()> r;
};



class NormalDist {

	public:       
        NormalDist(double mean, double stddev, double min, double max) : distribution(mean, stddev), min(min), max(max), r(std::bind(distribution, std::default_random_engine(seed))) {}
        double Generate() { double randNum = r(); 
        					if (randNum < min)
        						return min;
        					if (randNum > max)
        						return max; 
        					return randNum; }

    private:

    	double min;
    	double max;
        std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::normal_distribution<double> distribution;
        std::function<double()> r;
};



class ExponentialDist {
    public:

        ExponentialDist(double lambda, double min, double max) : distribution(lambda), min(min), max(max), r(std::bind(distribution, std::default_random_engine(seed))) {}
        double Generate() { double randNum = r(); 
        					if (randNum < min)
        						return min;
        					if (randNum > max)
        						return max; 
        					return randNum; }

    private:

    	double min;
    	double max;
        std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::exponential_distribution<double> distribution;
        std::function<double()> r;
};



class ChiSquaredDist {
    public:

        ChiSquaredDist(double degrees, double min, double max) : distribution(degrees), min(min), max(max), r(std::bind(distribution, std::default_random_engine(seed))) {}
        double Generate() { double randNum = r(); 
        					if (randNum < min)
        						return min;
        					if (randNum > max)
        						return max; 
        					return randNum; }

    private:

    	double min;
    	double max;
        std::mt19937::result_type seed = std::chrono::high_resolution_clock::now().time_since_epoch().count();
        std::chi_squared_distribution<double> distribution;
        std::function<double()> r;
};



#endif // RAND_H
