#ifndef TUMOR_GENERATOR_H
#define TUMOR_GENERATOR_H

#include <Eigen/Dense>
#include <random>
#include <iostream>
#include <thread>

class TumorGenerator { // a mutational network
    public:
    TumorGenerator(const Eigen::ArrayXd& param_initial_population,
                    const Eigen::ArrayXXd& param_transition_rates,
                    const double& param_tmax, const double& param_dt);

    // printing methods
    void print_parameters();
    void print_status();
    void print_thread();

    // tumor evolution methods
    void initiation();
    void evolve_step(); //One step forward

    // collect the data of a single realization
    // return a table of data points (a Eigen array)
    Eigen::ArrayXXd single_tumor(unsigned seed); // overkill fixed

    private:
    // tumor parameters: Markov transition matrix, initial condition, and end time 
    int ntype;
    Eigen::ArrayXd initial_population;
    Eigen::ArrayXXd transition_rates;
    double tmax;
    double dt;
   
    // tumor status: current population and current time
    Eigen::ArrayXd population;
    double time;

    // tumor randomness
    std::mt19937_64 rnd_generator;

    // random events
    double lifespan(const Eigen::ArrayXXd& weights);
    std::tuple<int,int> increment_type(const Eigen::ArrayXXd& weights);
};

#endif