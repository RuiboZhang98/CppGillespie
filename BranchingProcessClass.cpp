#include <iostream>
#include <iomanip>
#include <random>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <tuple>
#include <iterator> //for std::back_inserter

#include <future>
#include <vector>
#include <numeric> // for std::iota. std::reduce
#include <ranges> // for std::ranges
#include <execution> // for std::execution::par

#include <fstream>
#include <string>
#include "FileToMatrix.cpp"
#include <utility>

class TumorGenerator { // a mutational network
    public:
    template <typename DerivedA, typename DerivedB>
    TumorGenerator(const Eigen::ArrayBase<DerivedA>& param_initial_population,
                    const Eigen::ArrayBase<DerivedB>& param_transition_rates,
                    const double& param_tmax, const double& param_tgrid);

    // printing methods
    void print_parameters();
    void print_status();

    // tumor evolution methods
    void initiation(){
        time = 0;
        population = initial_population;
    }

    void evolve_step(); //One step forward
    double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights);
    template <typename Derived>
    std::tuple<int,int> increment_type(const Eigen::ArrayBase<Derived> &weights);

    // collect the data of a single realization
    // return a table of data points (a Eigen array)
    Eigen::ArrayXXd single_tumor(unsigned seed); // overkill fixed

    private:
    // tumor parameters: Markov transition matrix, initial condition, and end time 
    int ntype;
    Eigen::ArrayXd initial_population;
    Eigen::ArrayXXd transition_rates;
    double tmax;
    double tgrid;
   
    // tumor status: current population and current time
    Eigen::ArrayXd population;
    double time;

    // tumor randomness
    std::mt19937_64 rnd_generator;
};

template <typename DerivedA, typename DerivedB>
TumorGenerator::TumorGenerator(const Eigen::ArrayBase<DerivedA>& param_initial_population,
                    const Eigen::ArrayBase<DerivedB>& param_transition_rates,
                    const double& param_tmax, const double& param_tgrid){
    
    // collect parameters
    ntype = param_initial_population.size();
    initial_population = param_initial_population;
    transition_rates = param_transition_rates;
    tmax = param_tmax;
    tgrid = param_tgrid;
}

void TumorGenerator::print_status(){
    std::cout << " time = " << time << std::endl;
    std::cout << "current population = " << population.transpose() << std::endl;
}

void TumorGenerator::print_parameters(){
    // print tumor parameters
    std::cout << "initial population = " << initial_population.transpose() << std::endl; 
    std::cout << "transition rate matrix: \n" << transition_rates << std::endl;
    std::cout << "maximum time = " << tmax << std::endl;
    std::cout << "grid size of time (delta t) = " << tgrid << std::endl;
}

double TumorGenerator::lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights){
    // The function computes the waiting time for next population change.
    double lambda = weights.sum();
    std::exponential_distribution<> exp{ lambda };
    return exp(rnd_generator);
} 

template <typename Derived> //std::array<int,2>, don't need REF
std::tuple<int,int> TumorGenerator::increment_type(const Eigen::ArrayBase<Derived> &weights)
{
    // The function gives the type (a integer) that increases its population during 
    // population change.
    int ntype = weights.rows();
    std::discrete_distribution<> d(weights.reshaped().begin(),weights.reshaped().end());

    int random_entry = d(rnd_generator);
    int c = random_entry / ntype;
    int r = random_entry - c * ntype;

    int decrement = -1;
    if(r != c) decrement = r;

    return {decrement, c};
}

void TumorGenerator::evolve_step(){

    Eigen::ArrayXXd weights = transition_rates.colwise() * 
        Eigen::Map<Eigen::RowVectorXd>(population.data(), ntype).array().transpose();

    auto [decre_type, incre_type] = increment_type(weights); //edit 8/14

    population(incre_type) += 1;
    if (decre_type != -1) population(decre_type) -= 1;
    time += lifespan(weights);
}



// Get a single realization of a tumor 
Eigen::ArrayXXd TumorGenerator::single_tumor(unsigned seed){
    
    // set the seed of the random number generator
    rnd_generator.seed(seed);
    const int datalen = (int)(tmax / tgrid) + 1;  // length of recorded data
    const Eigen::VectorXd record_time = Eigen::VectorXd::LinSpaced(datalen, 0, tmax); // record population at these populations
    Eigen::ArrayXXd time_population = Eigen::ArrayXXd::Constant(datalen, ntype + 1, -1.0); // initiate the data table

    time_population.block(0, 0, datalen, 1) = record_time; // flush the first column with the record times

    // initialize a tumor
    initiation();
    Eigen::ArrayXd population_old = population;

    // data_index indicates the current row of the result
    int data_index = 0;

    while(data_index < datalen){
        for ( ; (data_index < datalen) && (record_time(data_index) <= time) ; data_index++){
                time_population.block(data_index, 1, 1, ntype) = population_old.transpose();
        }
        population_old = population;
        evolve_step();
    }

    
    return time_population;
}

// regular functions
// a function that reads parameters from files
std::tuple<double, double, int, Eigen::ArrayXd, Eigen::ArrayXXd> 
read_parameters(const string& txt_filepath, const string& csv_filepath)
{
    std::ifstream sim_param_file(txt_filepath);
    std::string name;
    double tmax, tgrid;
    int ntype, runs;

    sim_param_file >> name >> ntype;
    sim_param_file >> name >> tmax;
    sim_param_file >> name >> tgrid;
    sim_param_file >> name >> runs;

    FileToMatrix ftm(ntype, ntype + 1, csv_filepath);

    Eigen::ArrayXd initial_population = ftm().col(0).array();
    Eigen::ArrayXXd transition_rates(ftm().block(0, 1, ftm().rows(), ftm().cols() - 1));

    return std::make_tuple(std::move(tmax), std::move(tgrid), std::move(runs), 
                            std::move(initial_population), std::move(transition_rates));
}

// a single task that constructs a tumor and obtain a single realization
Eigen::ArrayXXd tumor_par(unsigned seed, double tmax, double tgrid, 
    Eigen::ArrayXd initial_population, Eigen::ArrayXXd transition_rates){

    // construct a tumor by the parameters
    TumorGenerator tumorObj1(initial_population, transition_rates, tmax, tgrid);
    Eigen::ArrayXXd data = tumorObj1.single_tumor(seed);

    return data;
}

int main()
{
    // In the main function, we use the task based concurrency to generate tumors.
    // The "tumors" could be treated as potiential outcomes of a single tumors, patients that suffer 
    //      from a same type of cancer with different symptoms, or tumors within a single person

    // read parameters 
    auto [tmax, tgrid, runs, initial_population, transition_rates] = read_parameters("ParamsSimulation.txt", "ParamsArray.csv");

    // assign seeds to realizations
    std::vector<unsigned> seeds(runs);
    std::iota(seeds.begin(), seeds.end(), 0);

    // store future objects
    std::vector<std::future<Eigen::ArrayXXd>> sample_tumor_par;
    for (auto k : seeds) sample_tumor_par.emplace_back(std::async(tumor_par, k, tmax, tgrid, initial_population, transition_rates));

    // store results (a vector of future Eigen arrays) //edit this part done 8/13
    std::vector<Eigen::ArrayXXd> population_result;
    population_result.reserve(sample_tumor_par.size());
    
    std::ranges::transform(sample_tumor_par, std::back_inserter(population_result), [](std::future<Eigen::ArrayXXd> &fut){
        return fut.get();
    });
    
    // store waiting time results
    std::vector<Eigen::ArrayXXd> waitingtime_result;
    waitingtime_result.reserve(population_result.size());

    for (auto array : population_result){
        waitingtime_result.push_back((array > 0).cast<double>());
    }


    // compute average population
    Eigen::ArrayXXd initial_array = Eigen::ArrayXXd::Zero(population_result[0].rows(), population_result[0].cols());
    
    Eigen::ArrayXXd average_population = 1.0 / runs * 
        std::reduce(std::execution::par, population_result.begin(), population_result.end(), initial_array);

    std::cout << "Average population:\n" << average_population << std::endl;
    
    // computing waiting time distribution
    Eigen::ArrayXXd waitingtime_dist = 1.0 / runs * 
        std::reduce(std::execution::par, waitingtime_result.begin(), waitingtime_result.end(), initial_array);


    waitingtime_dist.block(0, 0, waitingtime_dist.rows(), 1) = average_population.block(0, 0, average_population.rows(), 1);

    std::cout << "Waiting time distribution:\n" << waitingtime_dist << std::endl;

    return 0;
}