# CppGillespie
The main goal of this project is to simulate a multitype branching process. For more information about the mathematical task, please refer to [Sampling a multi-type branching process](GillespieAlg.md). I'm also giving a talk for this specific project on CppCon 2024. The slides are also going to be attached.

# Compile
To compile the cpp file, one first installs the [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page). Here the recommand compiler is the MSVC. Using this compiler, one can compile the code by the follwing command line (on a windows machine):

    cl -I /path/to/eigen/ /std:c++20 /Febranching.exe BranchingProcessPar.cpp .\TumorGenerator.cpp .\FileToMatrix.cpp

Here, 'branching' is the name of the executable file. 

One can also compile the code by the gcc compiler the following command line:

    g++ -I /path/to/eigen/ -std=c++20 BranchingProcessPar.cpp .\TumorGenerator.cpp .\FileToMatrix.cpp -o branching 

The code can be compiled, but there is going to be a warning saying: [Parallel STL message]: "Vectorized algorithm unimplemented, redirected to serial". This is because that the paralleled version of STL algorithm is currently (2024.9) only implemented in the MSVC compiler.

# Run the exectuable file

Use the following command to run the executable file:

    ./branching.exe

## Parameter values

The program is going to get the parameters by reading two files including ParamsSimulation.txt and ParamsArray.csv.

The ParasSimulation.txt is organized as

CML_Example ----------- output file name
ntype 3     ----------- the number of types
tmax 85.0   ----------- the maximum simulation time
tgrid 1     ----------- record the population of cells every tgrid
runs 1000   ----------- number of simulations

The ParamsArray.csv file contains a single matrix, e.g.

100000, 0, 1e-06, 0,
0, 0, 0.046, 1.1e-06,
0, 0, 0, 0,

Here the first column (100000,0,0) is the initial population. The rest of the matrix presents the transition rate matrix between all the types. More specifically, the diagonal entries are the birth rates. A off-diagonal entry $$A_ij$$ represents the mutation rate from type i to type j. 

# Output files

All the output files are collected in a single folder under the name CML_Example.

i. CML_Example_avgpopulation_runs1000tmax85

This file collects the average cellular population in the entire simulation time interval with a resolution $ \Delta t = tgrid$. It contains a 

$$\left(\frac{\text{tmax}}{\Delta t} + 1\right) \text{ by } \text{ number of types } + 1$$

matrix. The first column records times.

ii. CML_Example_params.txt

This files collects the initial condition and the transition rate matrix.

iii. CML_Example_waitingtime_runs1000tmax85

This files presents the waiting time distribution of all the types include the first type. The first column records times. It contains a

$$\left(\frac{\text{tmax}}{\Delta t} + 1\right) \text{ by } (\text{number of types } + 1)$$

matrix.
   

