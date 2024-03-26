# CppGillespie
The goal of this project is to build up a Gillespie simulation from the beginning using C plus plus. Along the way, all the issues are good practise of dealing with Cpp, Github, Git Bash, department server...

The main goal of this project is to simulate a multitype branching process.

# Compile
To compile the cpp file, one first installs the [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page). After having the Eigen library, one can compile the code by the following command line:

    g++ -I /path/to/eigen/ -std=c++20 CppGillespie.cpp -o CppGillespie 

One needs to replace "/path/to/eigen/" by the directory of Eigen.

# Run the exectuable file

Use the following command to run the executable file:

    ./CppGilespie -o CppGilespie --seed 1 --runs 1 --tmax 5

## Arguments

-o CppGilespie (string): This argument names output file names (Required)
--seed 1 (int): This argument provides the random seed in the random number generator (Required)
--runs 1 (positive int): This argument provides the number of independent samples. (Optional with default value 1)
--tmax 5 (positive double): This argument provides the maximum time in each realization. (Optional with default value 20.0)

## Special Option

Use `-debug` to print log file. To save the log in to a txt file, use

    ./CppGilespie -o CppGilespie --seed 1 --runs 1 --tmax 5 -debug >> log.txt

Note one needs to ==hit enter== after run the above line to write the log.

## Multithread Computaiton via bash script

Multithreadbash.sh is a simple bash script that allows the executable file to run parallelly on a multicore machine. 

    for sd in {21..40}
    do
        ./CppGilespie -o CppGilespie --seed $sd --runs 100 &
    done

In the above script, `for sd in {21..40}` is a for loop that assign different random seeds to different threads.



