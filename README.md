# impact-of-hypocrisy
Source code for simulations presented in "The impact of hypocrisy on opinion formation: a dynamic model"

## Building the executable
1. Install GCC (Linux/Windows) or XCode (Mac). 
2. Install the GNU Scientific Library, see https://www.gnu.org/software/gsl/.
3. In a terminal, run these commands.
```
git clone https://github.com/mgastner/impact-of-hypocrisy.git
cd impact-of-hypocrisy/
make
```
This will build two executables: `bvm` for the Basic Voter Model and `cvm` for the Concealed Voter Model.

## Using the executables

### Basic voter model
```
./bvm n_run rc n nr_init (seed)
```
where
- `n_run`: number of runs
- `rc`: rate of copying a neighbour
- `n`: number of agents
- `nr_init`: initial number of red agents
- `seed`: seed of the random number generator, automatically chosen as system time if left empty

### Concealed voter model
```
./cvm n_run rc re ri n nr_ext_init nr_int_init nrr_init (seed)
```
where
- `n_run`: number of runs
- `rc`: rate of copying a neighbour
- `re`: rate of externalization
- `ri`: rate of internalization
- `n`: number of agents
- `nr_ext_init`:  initial number of agents whose external opinion is red
- `nr_int_init`:  initial number of agents whose internal opinion is red
- `nrr_init`: initial number of agents whose opinion is red, both externally and internally.

## Output

- observed fraction *F* of runs in which the red opinion won
- predicted fraction *m* of runs in which the red opinion won
- observed mean consensus time with 95% confidence interval
- predicted mean consensus time
- During the first run, the fraction of red agents as a function of time is printed to the file `dynamics.dat`.
