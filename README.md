# CppGillespie
The goal of this project is to build up a Gillespie simulation from the beginning using C plus plus. Along the way, all the issues are good practise of dealing with Cpp, Github, Git Bash, department server...

The main goal of this project is to simulate a multitype branching process. For more information about the mathematical task, please refer to [Sampling a multi-type branching process](GillespieAlg.md).

# Compile
To compile the cpp file, one first installs the [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page). After having the Eigen library, one can compile the code by the following command line:

    g++ -I /path/to/eigen/ -std=c++20 CppGillespie.cpp -o CppGillespie 

One needs to replace "/path/to/eigen/" by the directory of Eigen.

# Run the exectuable file

Use the following command to run the executable file:

    ./CppGilespie -o CppGilespie --seed 1 --runs 5 --tmax 10

## Arguments

-o CppGilespie (string): This argument names output file names (Required)
--seed 1 (int): This argument provides the random seed in the random number generator (Required)
--runs 1 (positive int): This argument provides the number of independent samples. (Optional with default value 1)
--tmax 5 (positive double): This argument provides the maximum time in each realization. (Optional with default value 20.0)

## Options

1. Use `-debug` to print log file. To save the log in to a txt file, use

    ./CppGilespie -o CppGilespie --seed 1 --runs 5 --tmax 10 -debug >> log.txt

Note that one needs to **hit enter** after run the above line to write the log.

2. Use `-printpop` to print cell populations for all realizations. For example,

    ./CppGilespie -o CppGilespie --seed 1 --runs 5 --tmax 10 -printpop

is going to output an extra .txt file. Please refer to the description in the section below titled 'Output files'.

# Output files

i. (Optional) CppGilespie_population_seed1runs5tmax10.txt

This file collects the cellular population in the entire simulation time interval with a resolution $ \Delta t = 0.1$. It contains a 

$$\left(\text{runs} \times \left(\frac{\text{tmax}}{\Delta t} + 1\right)\right) \text{ by } \text{ number of types } + 1$$

matrix. The first column records times.

ii. CppGilespie_tau_seed1runs5.txt

This files collects the first arrival time of each type except the intial type in all the realizations. It contains a 

$$\text{runs} \text{ by } (\text{number of types } - 1)$$

matrix.

iii. CppGilespie_waitingtime_seed1runs5.txt

This files presents the waiting time distribution of all the types include the first type. The first column records times. It contains a

$$\text{runs} \text{ by } (\text{number of types } + 1)$$

matrix.

# Multithread computaiton via bash script

Multithreadbash.sh is a simple bash script that allows the executable file to run parallelly on a multicore machine. 

    for sd in {21..40}
    do
        ./CppGilespie -o CppGilespie --seed $sd --runs 100 &
    done

In the above script, `for sd in {21..40}` is a for loop that assign different random seeds to different threads.

# Some future plans

1. Add OOP stuff
2. Read file for input model parameters
3. Figure out the correct way to do the parallel computing inside c++ (instead of a bash script)
   

