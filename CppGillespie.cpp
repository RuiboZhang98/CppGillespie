#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

#include <Eigen/Dense> // path: G:\\code\\cfiles\\eigen-3.4.0
#include <Eigen/Core>  // C:\\Users\\rzhang98\\source\\library\\eigen-3.4.0

using std::vector, std::string;
using Eigen::Ref, Eigen::ArrayXXd, Eigen::ArrayXd, Eigen::ArrayXi;

// function prototypes
double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng);

template <typename Derived>
int increment_type(const Eigen::ArrayBase<Derived> &weights, std::mt19937_64 &rng, int &decrement);



int main(int argc, char** argv)
{
    using std::cout, std::endl;
    using Eigen::Map, Eigen::RowVectorXd, Eigen::VectorXd, Eigen::seq; //, Eigen::last;

    // default parameters that are related to command line arguments
    string outputfile = "";             // related to the output file name
    long long unsigned int seed = 0;    // random seed
    int runs = 1;                       // number of realizations
    double tmax = 20;                   // end time of simulations
    bool debug_mode = false;            // debug mode prints traj
    bool printpop = false;

    // read parameters
    for (int i = 1; i < argc; ++i) {
        // Errors
        if ((string)argv[i] == "-o") {
            if (i == argc-1) {
                cout << "Error: no output file" << endl;
                return 1;
            } else {
                // set the name (or common theme) of the output file
                outputfile = (string)argv[i+1];
            }
        }

        if ((string)argv[i] == "--seed") {
            if (i == argc-1) {
                cout << "Error: please enter the seed for generating random numbers" << endl;
                return 1;
            } else {
                // set the random number seed
                seed = (long long unsigned int)atoi(argv[i+1]);
            }
        }

        if ((string)argv[i] == "--runs") {
            if (i == argc-1) {
                cout << "Error: please enter the number of runs" << endl;
                return 1;
            } else {
                // set number of runs
                runs = atoi(argv[i+1]);
            }
        }

        if ((string)argv[i] == "--tmax") {
            if (i == argc-1) {
                cout << "Error: please enter the end time" << endl;
                return 1;
            } else {
                // set end time, default value is 20
                tmax = atoi(argv[i+1]);
            }
        }

        if ((string)argv[i] == "-debug") {
            // if I want to print the log
            debug_mode = true;
        }

        if ((string)argv[i] == "-printpop") {
            printpop = true;
        }

        if (outputfile == "") {
            cout << "Error: no output file" << endl;
            return 1;
        }
    }

    // branching process parameters
    const ArrayXd initial_population { {11000, 0, 0} }; 
    const int ntype = initial_population.size();
    const ArrayXd birth_rates { {0, 0.06, 0} };
    // ArrayXd death_rates { {0, 0, 0, 0 } }; We don't consider death now.
    const double u = 1e-6;          // The scale of mutation rates
    ArrayXXd transition_rates {      // This gives the ratio of  
        {0, 2, 0},             // the mutation rates (from type i to type j).
        {0, 0, 2.5},               // We need more computations to reach the actual mutation rates
        {0, 0, 0}
    };

    transition_rates = transition_rates * u;

    for (int i = 0; i < ntype; ++i){
        transition_rates(i, i) = birth_rates(i);
    }

    cout << "number of types = " << ntype << endl; // print some important parameters
    cout << "transition rate matrix = \n" << transition_rates << endl;


    // Simulation Parameters
    const double tgrid =  0.1;    // time grid length. On grid points populations are recorded
    const int datalen = (int)(tmax / tgrid) + 1;  // length of recorded data

    // output file names
    string outputfile_population_filename = outputfile + "_population_seed" 
    + std::to_string((int)seed) + "runs" + std::to_string(runs) + "tmax" + 
    std::to_string((int)tmax) + ".txt";
    
    string outputfile_waitingtime_filename = outputfile + "_waitingtime_seed" 
    + std::to_string((int)seed) + "runs" + std::to_string(runs) + ".txt";

    string outputfile_tau_filename = outputfile + "_tau_seed" 
    + std::to_string((int)seed) + "runs" + std::to_string(runs) + ".txt";


    cout << "number of runs = " << runs << endl; // print some important parameters
    cout << "tmax = " << tmax << endl;
    cout << "tgrid = " << tgrid << endl;
    cout << "printpop = " << printpop << endl;

    std::mt19937_64 mt { seed };        // random number generator. When having concurrency, the seed for each thread should be different 
    
    const VectorXd record_time = VectorXd::LinSpaced(datalen, 0, tmax);
    ArrayXd population(ntype);
    ArrayXd old_population(ntype); // record the population of the previous time step
    ArrayXXd weights(ntype, ntype);

    // In the following population data array, each colume stores a single realization.
    // The first column is used to store the time grid;
    ArrayXXd population_data = ArrayXXd::Constant(datalen * runs, ntype + 1,-1.0);
    
    // The following matrix is used to count if type i exists at record_time points
    ArrayXXd waitingtime_data = ArrayXXd::Constant(datalen, ntype + 1, 0.0);

    // The following matrix is used to collect the actual arrival time
    ArrayXXd tau_data = ArrayXXd::Constant(runs, ntype - 1, -1.0);    

    if(printpop)
        population_data.block(0, 0, datalen, 1) = record_time; // write times in the first column
    
    waitingtime_data.block(0, 0, datalen, 1) = record_time;

    int data_index, run_index, incre_type, decre_type;
    double t;

    if (debug_mode){
        cout << endl << "record_time: " << record_time << endl;
        cout << "Press enter to continue " << endl;
        getchar();
    }

    // The main body of the Gillespie Simulation
    for (run_index = 0; run_index < runs; ++run_index){
        population = initial_population;     // reset population
        old_population = initial_population;
        t = 0;              // reset time
        data_index = 0;     // reset the index of data

        while (t < tmax){
            for ( ; record_time(data_index) < t ; ++data_index){
                if(printpop){
                    population_data.block(run_index * datalen + data_index, 1, 1, ntype) = old_population.transpose();
                    population_data(run_index * datalen + data_index, 0) = record_time(data_index);
                }
                    
                waitingtime_data.block(data_index, 1, 1, ntype) += (old_population.transpose() > 0).cast<double>();
                waitingtime_data(data_index, 0) = record_time(data_index);

                if (debug_mode){
                    cout << "t = " << t 
                    << " record_time = " << record_time(data_index)
                    << " Population = " << population.transpose() << endl;
                }
            }

            weights = transition_rates.colwise() * Map<RowVectorXd>(population.data(), population.size()).array().transpose();
            
            decre_type = -1;
            incre_type = increment_type(weights, mt, decre_type);
            old_population = population;
            population(incre_type) += 1;
            if (decre_type != -1)
                population(decre_type) -= 1;

            t+= lifespan(weights, mt);

            if (incre_type != 0 && population(incre_type) == 1)
                tau_data(run_index, incre_type - 1) = t;

            if (debug_mode){
                cout << "t = " << t << "old: " << old_population.transpose() << 
                "currect: " << population.transpose() << endl;
            }
        }
        // for the end time
        for ( ; record_time(data_index) < t ; ++data_index){
            if(printpop){
                population_data.block(run_index * datalen + data_index, 1, 1, ntype) = old_population.transpose();
                population_data(run_index * datalen + data_index, 0) = record_time(data_index);
            }
            
            waitingtime_data.block(data_index, 1, 1, ntype) += (old_population.transpose() > 0).cast<double>();
            waitingtime_data(data_index, 0) = record_time(data_index);

            if (debug_mode){
                cout << "End Time: t = " << t
                << " record_time = " << record_time(data_index)
                << " Population = " << population.transpose() << endl;

                cout << "data_index = " << data_index << 
                "record size = " << record_time.size() << "datalen = " 
                << datalen << endl;
            }
            
            if (data_index + 1 == record_time.size())
                break;
        }

        cout << "The " << run_index + 1 << "th run finishes" << endl;
    }
    // normalize waiting time data

    waitingtime_data.block(0, 1, datalen, ntype) /= runs;

    // save to files
    std::ofstream outFile;

    if(printpop){
        outFile.open(outputfile_population_filename);
        outFile << endl << population_data << endl;
        outFile.close();
    }

    outFile.open(outputfile_waitingtime_filename);
    outFile << endl << waitingtime_data << endl;
    outFile.close();

    outFile.open(outputfile_tau_filename);
    outFile << endl << tau_data << endl;
    outFile.close();

    return 0;
}

// function definitions
double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng)
{
    // The function computes the waiting time for next population change.
    double lambda = weights.sum();
    std::exponential_distribution<> exp{ lambda };
    return exp(rng);
} 

template <typename Derived>
int increment_type(const Eigen::ArrayBase<Derived> &weights, std::mt19937_64 &rng, int &decrement)
{
    // The function gives the type (a integer) that increases in a population change.
    int ntype = weights.rows();
    std::discrete_distribution<> d(weights.reshaped().begin(),weights.reshaped().end());

    int random_entry = d(rng);
    int c = random_entry / ntype;
    int r = random_entry - c * ntype;

    if(r != c) decrement = r;

    return c;
}

