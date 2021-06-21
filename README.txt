The DOME modeling tool is implemented in MATLAB. It requires the Matlab Control System toolbox and the CVX toolbox.

Lin_sys = dome (chunk,Ts,r,prec,epsilon) returns the transfer function representation of the linear system identified by the algorithm. The structure of inputs to the "dome" function is as follows: 
1. "chunk" is a structured array whose size is the number of chunks and where each of the entries has two elements in its structure, "chunk(i).stim" and "chunk(i).out". These are two arrays of the same size with the chunk(i).stim containing the values of the
2. "Ts" is the sampling period of the data, which is equal for all chunks.
3. "r" is the degree of discretization of the unit circle and discerns the total number of candidate poles that is approximately proportional to r^2.
4. "prec" is a bound used to decide when a quantity is assumed to be zero. In other
5. "epsilon" is the small constant to make the optimization problem well-posed.

Two datasets named "data_set_1" and "data_set_2" are available. Loading each dataset in MATLAB loads two variables into the workspace: "chunk" and "Ts".