function [chunk] = chunk_creator(X,Y,len_chunks,ll)

% This function creates structure arrays of stimulus and ooutcome.
% The inputs to this function are X that is the stimulus samples,
% Y that is the outcome samples, len_chunks that is a vector containing
% all number of samples in each chunk, and ll that is the total number of
% chunks.
% The output of the function, chunk, is a data structure containig the
% stimulus/outcome samples structured in chunks, where "stim" has the
% stimulus data and "out" contains the outcome data.


i = 1;
j = len_chunks(1);
len_chunks(1+ll) = 0;
for l=1:ll
    stim= X(i:j,1);
    out = Y(i:j,1);
    chunk(l).stim = stim;
    chunk(l).out = out;
    i = j+1;
    j = j+ len_chunks(l+1);
end