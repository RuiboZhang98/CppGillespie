# CppGillespie
The goal of this project is to build up a Gillespie simulation from the beginning using C plus plus. Along the way, all the issues are good practise of dealing with Cpp, Github, Git Bash, department server...

The main goal of this project is to simulate a multitype branching process.

# Compile
To compile the cpp file, one first install the [Eigen library](https://eigen.tuxfamily.org/index.php?title=Main_Page). After having the Eigen library, one can compile the code by the following command line:

    g++ -I /path/to/eigen/ CppGillespie.cpp -o CppGillespie

One needs to replace "/path/to/eigen/" by the directory of Eigen.

# Run the .exe
Simply use ./CppGillespie to run the code. It will create a file named "data.txt" in the current directory. 
