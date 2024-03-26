#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

#include <Eigen/Dense> // path: G:\\code\\cfiles\\eigen-3.4.0

using std::vector, std::string;
using Eigen::Ref, Eigen::ArrayXXd, Eigen::ArrayXd, Eigen::ArrayXi;

// function prototypes
int find_positive_indices(vector<int>* r, vector<int>* c, 
const Eigen::Ref<const Eigen::ArrayXXd>& matrix, const int& ntype);
double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng);
int increment_type(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng, int &flag);

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
                // set number of threads
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
            debug_mode = true;
        }

        if (outputfile == "") {
            cout << "Error: no output file" << endl;
            return 1;
        }
    }

    // branching process parameters
    const ArrayXd initial_population { {10000, 0, 0, 0, 0} }; 
    const int ntype = initial_population.size();
    const ArrayXd birth_rates { {0, 0, 0.7, 1.0, 1.0} };
    // ArrayXd death_rates { {0, 0, 0, 0 } }; We don't consider death now.
    const double u = 1e-3;          // The scale of mutation rates
    ArrayXXd transition_rates {      // This gives the ratio of  
        {0, 0.2, 0, 0, 0},             // the mutation rates (from type i to type j).
        {0, 0, 8, 0, 0},               // We need more computations to reach the actual mutation rates
        {0, 0, 0, 5, 0},
        {0, 0, 0, 0, 6},
        {0, 0, 0, 0, 0}
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
    int flag;

    // output file name
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

    std::mt19937_64 mt { seed };        // random number generator. When having concurrency, the seed for each thread should be different 
    
    const VectorXd record_time = VectorXd::LinSpaced(datalen, 0, tmax);
    ArrayXd population(ntype);
    ArrayXd old_population(ntype); // record the population of the previous time step
    ArrayXXd weights(ntype, ntype);

    // In the following population data array, each colume stores a single realization.
    // The first column is used to store the time grid;
    ArrayXXd population_data = ArrayXXd::Constant(datalen * runs, ntype + 1,-1.0);
    
    // The following matrix is used to count if type i exists at record_time points
    ArrayXXd waitingtime_data = ArrayXXd::Constant(datalen, ntype + 1, -1.0);

    // The following matrix is used to collect the actual arrival time
    ArrayXXd tau_data = ArrayXXd::Constant(runs, ntype - 1, -1.0);    

    population_data.block(0, 0, datalen, 1) = record_time; // write times in the first column
    waitingtime_data.block(0, 0, datalen, 1) = record_time;

    int data_index, run_index, change_index;
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
                population_data.block(run_index * datalen + data_index, 1, 1, ntype) = old_population.transpose();
                waitingtime_data.block(data_index, 1, 1, ntype) += (old_population.transpose() > 0).cast<double>();

                population_data(run_index * datalen + data_index, 0) = record_time(data_index);
                waitingtime_data(data_index, 0) = record_time(data_index);

                if (debug_mode){
                    cout << "t = " << t 
                    << " record_time = " << record_time(data_index)
                    << " Population = " << population.transpose() << endl;
                }
            }

            weights = transition_rates.colwise() * Map<RowVectorXd>(population.data(), population.size()).array().transpose();
            change_index = increment_type(weights, mt, flag);
            old_population = population;
            population(change_index) += 1;
            if (flag == 1)
                population(change_index - 1) -= 1;
            t+= lifespan(weights, mt);

            if (change_index != 0 && population(change_index) == 1)
                tau_data(run_index, change_index - 1) = t;

            if (debug_mode){
                cout << "t = " << t << "old: " << old_population.transpose() << 
                "currect: " << population.transpose() << endl;
            }
        }
        // for the end time
        for ( ; record_time(data_index) < t ; ++data_index){
            population_data.block(run_index * datalen + data_index, 1, 1, ntype) = old_population.transpose();
            waitingtime_data.block(data_index, 1, 1, ntype) += (old_population.transpose() > 0).cast<double>();
            population_data(run_index * datalen + data_index, 0) = record_time(data_index);
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
    waitingtime_data /= runs;

    // save to files
    std::ofstream outFile;
    outFile.open(outputfile_population_filename);
    outFile << endl << population_data << endl;
    outFile.close();

    outFile.open(outputfile_waitingtime_filename);
    outFile << endl << waitingtime_data << endl;
    outFile.close();

    outFile.open(outputfile_tau_filename);
    outFile << endl << tau_data << endl;
    outFile.close();

    return 0;
}

// function definitions

int find_positive_indices(vector<int>* r, vector<int>* c, const Eigen::Ref<const Eigen::ArrayXXd>& matrix, const int& ntype)
{
    // This function finds the positive indices of a given 2-d Eigen Array.
    // The row indices and column indicies are recorded in vectors r and c, respectively.
    // The function returns the number of positive elements.
    int positive_num = (matrix > 0).count();
    if (positive_num > 0){
        for (int i = 0; i < ntype ; ++i){
            for (int j = 0; j < ntype ; ++j){
                if (matrix(i,j) > 0){
                    // std::cout << i << ", " << j << ", " << matrix(i,j) << std::endl;
                    r->push_back(i);
                    c->push_back(j);
                }
            }
        }
    }
    return positive_num;
}

double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng)
{
    // The function computes the waiting time for next population change.
    double lambda = weights.sum();
    std::exponential_distribution<> exp{ lambda };
    return exp(rng);
} 

int increment_type(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng, int &flag)
{
    // The function gives the type (a integer) that increases in a population change.
    int ntype = weights.rows();
    std::uniform_real_distribution<> unif { 0.0, 1.0 };

    vector<int> r, c, positive_entries;
    int positive_num = find_positive_indices(&r, &c, weights, ntype); // find possible events
    ArrayXXd prob = weights / weights.sum();
    // std::cout << prob << std::endl;
    double single_draw = unif(rng);
    int i;
    for(i = 0; i < positive_num; ++i){
        single_draw -= prob(r[i],c[i]);
        if (single_draw < 0)
            break;
    }
    flag = 1;
    if (r[i] == c[i])
        flag = 0;
    return c[i];
}