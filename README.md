# CX-4220-PA-1
Programming Assignment 1

Setting up git:

Download git via command line tools

$ git config --global user.name "YOUR NAME"

$ git config --global user.email "YOUR EMAIL USED FOR GITHUB"

Navigate to where you want to clone the repository (where you want the folder containing program files to go).
This part might be weird

$ git clone https://github.gatech.edu/ttippens6/CX-4220-PA-1.git

You now have a local copy you can work on.

Using git:

Before making changes, do

$ git pull

$ git checkout master

Then do your changes.
Once your changes have been added, you 'commit them'. Commit once, after your done working.
Do:

$ git commit -a -m "Brief description of the changes you made"

After the commit, you'll need to push your changes to the online repo:

$ git push

and you're done.

HPC [Spring 17] - Programming Assignment 1
========================================

## Code hosting

You are highly encouraged to use a code version management tool such as git.
This will help you to code in a team and keep track of your progress.

However, do not upload your code to public repositories. If you want to use
version management for collaboration, make sure to use private repositories.

We highly recommend using Georgia Tech's Enterprise Github installation at
https://github.gatech.edu/

We will also be hosting the framework code for the programming assignment on
there.  If you find any issues with the code framework and we have to make
changes, we will publish those changes in that GitHub repository additionally to
sending out the updated framework.

## Code structure

All the code is located at the root level of the project.

There are multiple header and .cpp files, your implementation will go
into the following files:

- `evaluator.cpp`: Implement the sequential algorithm for polynomial evaluation
  to the function declaration in `evaluator.h`
- `mpi_evaluator.cpp`: Implement the parallel algorithm for polynomial evaluation
  to the function declaration in `mpi_evaluator.h`
- `utils.cpp`: Implement any utility functions you need while implementing your
  algorithms put the function declarations in `utils.h`

Other files containing code that you should not change are:

- `main.cpp`: Implements code for the main executable `poly-eval`. This does
  input/output reading and calling of the actual functions.
- `io.h`: implements IO functions and random input generation
- `const.h`: contains constants used in the application

## Compiling

In order to compile everything, simply run
```sh
make all
```
## Running
### Serial execution example:
```sh
./poly-eval sample-constants.txt sample-values.txt
```
### Parallel execution example:
```sh
mpirun -np 3 ./poly-eval sample-constants.txt sample-values.txt
```
