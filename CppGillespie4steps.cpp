#include <iostream>
#include <fstream>
#include <vector>
#include <random>

#include <Eigen/Dense>

using std::vector;
using Eigen::Ref, Eigen::ArrayXXd, Eigen::ArrayXd, Eigen::ArrayXi;

// function prototypes
int find_positive_indices(vector<int> * r, vector<int> * c, 
const Eigen::Ref<const Eigen::ArrayXXd>& matrix, const int& ntype);

double lifespan(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng);
int increment_type(const Eigen::Ref<const Eigen::ArrayXXd>& weights, std::mt19937_64& rng, int &flag);

int main()
{
    using std::cout, std::endl;
    using Eigen::Map, Eigen::RowVectorXd, Eigen::VectorXd, Eigen::seq; //, Eigen::last;

    // branching process parameters
    const ArrayXd initial_population { {10000, 0, 0, 0 } }; 
    const int ntype = initial_population.size();
    const ArrayXd birth_rates { {0, 0, 1, 1.5 } };
    // ArrayXd death_rates { {0, 0, 0, 0 } }; We don't consider death now.
    const double u = 1e-3;          // The scale of mutation rates
    ArrayXXd transition_rates {      // This gives the ratio of  
        {0, 2, 0, 0},             // the mutation rates (from type i to type j).
        {0, 0, 1, 0},               // We need more computations to reach the actual mutation rates
        {0, 0, 0, 5},
        {0, 0, 0, 0}
    };

    transition_rates = transition_rates * u;

    for (int i = 0; i < ntype; ++i){
        transition_rates(i, i) = birth_rates(i);
    }

    cout << "number of types = " << ntype << endl; // print some important parameters
    cout << "transition rate matrix = \n" << transition_rates << endl;

    // Simulation Parameters
    const double tmax = 20;          // ending time of simulations
    // have to add a little bit time to actually reach tmax
    const int runs = 1;        // number of realizations
    const double tgrid =  0.1;    // time grid length. On grid points populations are recorded
    const int datalen = (int)(tmax / tgrid);  // length of recorded data
    int flag;

    cout << "number of runs = " << runs << endl; // print some important parameters
    cout << "datalen = " << datalen << endl;

    long long unsigned int seed = 0;    // random seed
    std::mt19937_64 mt { seed };        // random number generator. When having concurrency, the seed for each thread should be different 
    
    const VectorXd record_time = VectorXd::LinSpaced(datalen, 0, tmax);
    ArrayXd population(ntype);
    ArrayXXd weights(ntype, ntype);
    // In the following population data array, each colume stores a single realization.
    // The first column is used to store the time grid;
    ArrayXXd population_data = ArrayXXd::Constant(datalen * runs, ntype + 1,-1.0); 
    
    population_data.block(0, 0, datalen, 1) = record_time; // write times in the first column

    int data_index, run_index, change_index;
    double t;
    // save to a file
    std::ofstream outFile;
    outFile.open("data.txt");

    // The main body of the Gillespie Simulation
    for (run_index = 0; run_index < runs; ++run_index){
        population = initial_population;     // reset population
        t = 0;       // reset time
        data_index = 0;     // reset the index of data
        while (t < tmax){
            for ( ; record_time(data_index) <= t ; ++data_index){
                population_data.block(run_index * datalen + data_index, 1, 1, 4) = population.transpose();
                cout << "record_time = " << record_time(data_index)
                << " Population = " << population.transpose() << endl;
            }
            weights = transition_rates.colwise() * Map<RowVectorXd>(population.data(), population.size()).array().transpose();
            change_index = increment_type(weights, mt, flag);
            population(change_index) += 1;
            if (flag == 1)
                population(change_index - 1) -= 1;
            t+= lifespan(weights, mt);
            // cout << "t = " << t << endl;
        }
        // for the end time
        cout << "record_time = " << record_time(data_index)
            << " Population = " << population.transpose() << endl;
        population_data.block(run_index * datalen + data_index, 1, 1, 4) = population.transpose();
        cout << "The " << run_index + 1 << "th run finishes" << endl;
    }

    outFile << endl << population_data << endl;
    outFile.close();

    cout << "Finished. Press any key to exit.";
    getchar();

    return 0;
}

// function definitions

int find_positive_indices(vector<int> * r, vector<int> * c, const Eigen::Ref<const Eigen::ArrayXXd>& matrix, const int& ntype)
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