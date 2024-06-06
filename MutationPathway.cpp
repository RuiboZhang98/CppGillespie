#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <vector>

#include "GeneClass.h"

class MutationPathway { // a linear mutational pathway
    public:
    MutationPathway(const std::vector<Gene>& genes, int initial_population);
    void print_mutation_path_info();

    protected: 
    int ntype;
    std::vector<Gene> genes;
    Eigen::ArrayXd initial_population, birth_rates;
    Eigen::ArrayXXd transition_rates;
};

class Simulator: public MutationPathway { // Gillespie Simulation
    public:
    Simulator(int runs, int seed, float tmax, float tgrid);
    void print_simulator_info();

    private: 
    // simulation parameters
    int runs, seed;
    float tmax;

    // when an how to record data
    float tgrid;
    int datalen;
    bool printpop;
    Eigen::VectorXd record_time;

    // output file
    std::string filePathAndName;
};

int main(){

    using std::vector;
    using Eigen::Ref, Eigen::ArrayXXd, Eigen::ArrayXd, Eigen::RowVectorXd;

    TumorSuppressor APC("genes/G1TSGApc.txt"), TP53("genes/G3TSGTp53.txt");
    Ocogene Kras("genes/G2OcogeneKras.txt");
    
    APC.print_gene_info();
    Kras.print_gene_info();
    TP53.print_gene_info();

    vector<Gene> genes = {APC, Kras, TP53};

    return 0;
}