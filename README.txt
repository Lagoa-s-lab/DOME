This package is implemented in MATLAB. It requires the Matlab Control System toolbox and the CVX toolbox.

The dome.m file returns the identified system with the following specifications:
Lin_sys = dome (chunk,Ts,r,prec,epsilon) returns the transfer function representation of the linear system identified by the algorithm. The structure of inputs to the "dome" function is as follows: 
1. "chunk" is a structured array whose size is the number of chunks and where each of the entries has two elements in its structure, "chunk(i).stim" and "chunk(i).out". These are two arrays of the same size with the chunk(i).stim containing the values of the
stimulus for the i-th chunk and chunk(i).out containing the corresponding outcomes. Each entry of chunk(i).out takes the value 1 if the outcome is true, 0 if it is not true and -1 if outcome is not available.
2. "Ts" is the sampling period of the data, which is equal for all chunks.
3. "r" is the degree of discretization of the unit circle and discerns the total number of candidate poles that is approximately proportional to r^2.
4. "prec" is a bound used to decide when a quantity is assumed to be zero. In other
words, if |c| <= prec then we assume that c = 0.
5. "epsilon" is the small constant to make the optimization problem well-posed.

The dome_train_test.m file is a similar function as dome but instead "chunk" as the input, it gets "chunk_train" and "chunk_test" and returns the error between the estimated and the true outcomes in addition to the identified system, Lin_sys. 

Four datasets, "data_set_1" and "data_set_2", "data_set_3", "data_set_4" are also provided with this package. Loading the first three datasets in MATLAB results in two variables into the workspace: "chunk" and "Ts" that are the first two inputs to "dome" function. data_set_1 is a small dataset with 100 samples and a simple model, data_set_2 is larger and has a more complex model. data_set_3 is a large set of 1999 samples with 50 chunks and 20% missing data in the outcome.
data_set_4 contains chunk_train, chunk_test, and Ts as the three first inputs to the function "dome_train_test". chunk_train consists of 45 chunks and chunk_train contains 5 chunks and the number of samples are not evenly distributed. Missing data is present in both training and test sets. 

dome_response_plots.m provides the code to plot impulse and pulse response (response to a 10-min pulse) of the system.
