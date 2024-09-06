#include <filesystem>

#include <future>
#include <vector>
#include <numeric> // for std::iota. std::reduce
#include <ranges> // for std::ranges
#include <execution> // for std::execution::par

#include "TumorGenerator.h"
#include "FileToMatrix.h"

// regular functions
// a function that reads parameters from files
std::tuple<std::string, double, double, int, Eigen::ArrayXd, Eigen::ArrayXXd> 
read_parameters(const std::string& txt_filepath, const std::string& csv_filepath)
{
    std::ifstream sim_param_file(txt_filepath);
    std::string name, filename;
    double tmax, dt;
    int ntype, runs;

    sim_param_file >> filename;
    sim_param_file >> name >> ntype;
    sim_param_file >> name >> tmax;
    sim_param_file >> name >> dt;
    sim_param_file >> name >> runs;

    FileToMatrix ftm(ntype, ntype + 1, csv_filepath);

    Eigen::ArrayXd initial_population = ftm().col(0).array();
    Eigen::ArrayXXd transition_rates(ftm().block(0, 1, ftm().rows(), ftm().cols() - 1));

    return std::make_tuple(std::move(filename), std::move(tmax), std::move(dt), std::move(runs), 
                            std::move(initial_population), std::move(transition_rates));
}

// a single task that constructs a tumor and obtain a single realization
Eigen::ArrayXXd tumor_par(const unsigned& seed, const double& tmax, const double& dt, 
    const Eigen::ArrayXd& initial_population, const Eigen::ArrayXXd& transition_rates){

    // construct a tumor by the parameters
    TumorGenerator tumorObj1(initial_population, transition_rates, tmax, dt);
    // tumorObj1.print_thread();
    Eigen::ArrayXXd data = tumorObj1.single_tumor(seed);
    // std::cout << "Seed " << seed << " Running on processor " << GetCurrentProcessorNumber() << std::endl;

    return data;
}

// write tables into txt files
void output_results(std::string filename, const double& tmax, const int& runs, 
                    const Eigen::ArrayXd& initial_population, 
                    const Eigen::ArrayXXd& transition_rates,
                    const Eigen::ArrayXXd& avg_population,
                    const Eigen::ArrayXXd& waiting_time_dist){
    
    
    std::filesystem::path path = filename;
    std::filesystem::create_directory(path);
    
    filename = filename + "\\" + filename;
    
    std::string waitingtime_filename = filename + "_waitingtime_runs" 
    + std::to_string(runs) + "tmax"  + std::to_string((int)tmax) + ".txt";

    std::string popuation_filename = filename + "_avgpopulation_runs" 
    + std::to_string(runs) + "tmax"  + std::to_string((int)tmax) + ".txt";

    std::string parameter_filename = filename + "_params.txt";
    
    std::ofstream outFile;
    outFile.open(parameter_filename);
    outFile << initial_population.transpose() << std::endl << transition_rates;
    outFile.close();

    outFile.open(waitingtime_filename);
    outFile << waiting_time_dist;
    outFile.close();

    outFile.open(popuation_filename);
    outFile << avg_population;
    outFile.close();
}

int main()
{
    // In the main function, we use the task based concurrency to generate tumors.
    // The "tumors" could be treated as potiential outcomes of a single tumors, patients that suffer 
    //      from a same type of cancer with different symptoms, or tumors within a single person

    // read parameters 
    auto [filename, tmax, dt, runs, initial_population, transition_rates] 
        = read_parameters("ParamsSimulation.txt", "ParamsArray.csv");

    // assign seeds to realizations
    std::vector<unsigned> seeds(runs);
    std::iota(seeds.begin(), seeds.end(), 0);

    // store future objects
    std::vector<std::future<Eigen::ArrayXXd>> sample_tumor_par;
    for (auto k : seeds) 
        sample_tumor_par.emplace_back(std::async(std::launch::async, tumor_par, k, tmax, dt, initial_population, transition_rates));

    // store results (a vector of future Eigen arrays) //edit this part done 8/13
    std::vector<Eigen::ArrayXXd> population_result;
    population_result.reserve(sample_tumor_par.size());
    
    std::ranges::transform(sample_tumor_par, std::back_inserter(population_result), [](std::future<Eigen::ArrayXXd> &fut){
        return fut.get();
    });

    // store waiting time results
    std::vector<Eigen::ArrayXXd> waitingtime_result;
    waitingtime_result.reserve(population_result.size());

    std::ranges::transform(population_result, std::back_inserter(waitingtime_result),
        [](const Eigen::ArrayXXd& array)->Eigen::ArrayXXd {return (array > 0).cast<double>();});

    // compute average population
    Eigen::ArrayXXd initial_array = Eigen::ArrayXXd::Zero(population_result[0].rows(), population_result[0].cols());
    
    Eigen::ArrayXXd average_population = 1.0 / runs * 
        std::reduce(std::execution::par, population_result.begin(), population_result.end(), initial_array);
    
    // computing waiting time distribution
    Eigen::ArrayXXd waitingtime_dist = 1.0 / runs * 
        std::reduce(std::execution::par, waitingtime_result.begin(), waitingtime_result.end(), initial_array);

    waitingtime_dist.block(0, 0, waitingtime_dist.rows(), 1) = average_population.block(0, 0, average_population.rows(), 1);

    output_results(filename, tmax, runs, initial_population, transition_rates, average_population, waitingtime_dist);

    return 0;
}