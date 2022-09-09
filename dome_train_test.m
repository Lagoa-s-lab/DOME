function [Lin_sys, error] = dome_train_test (chunk_test,chunk_train,Ts,r,prec,epsilon)

% This function takes 6 inputs, chunk_test,chunk_train,Ts,r,prec,epsilon, and returns 2 outputs, Lin_sys and error.
% [Lin_sys, error] = dome_train_test (chunk_test,chunk_train,Ts,r,prec,epsilon) returns the transfer function representation of 
% the linear system identified by the algorithm and the error between the estimated and the true outcomes.
% The structure of inputs to the "dome" function is as follows: 
% 1. "chunk_test" and "chunk_train" are structured arrays whose size are the number of chunks and where each of the entries has two elements in its structure,
% "chunk(i).stim" and "chunk(i).out". These are two arrays of the same size with the chunk(i).stim containing the values of the
% stimulus for the i-th chunk and chunk(i).out containing the corresponding outcomes in each data structure.
% Each entry of chunk(i).out takes the value 1 if the outcome is true, 0 if it is not true and -1 if outcome is not available.
% 2. "Ts" is the sampling period of the data, which is equal for all chunks.
% 3. "r" is the degree of discretization of the unit circle and discerns the total number of candidate poles that is approximately proportional to r^2.
% 4. "prec" is a bound used to decide when a quantity is assumed to be zero. In other
% words, if |c| <= prec then we assume that c = 0.
% 5. "epsilon" is the small constant to make the optimization problem well-posed.


% Creating poles consisting complex and real poles, with P being the number of poles
poles_real = linspace(-1,1,r);
poles_vec = [];
for rp = 1:r
    poles_vec = [poles_vec,poles_real(rp)+1j*linspace(0.1,1,r)];
end
poles_vec = [poles_vec,conj(poles_vec)];
poles_vec = poles_vec(find(abs(poles_vec) < 1));
P_complex = length(poles_vec);
poles_vec = [poles_vec,linspace(-1,1,2*r)];
poles_vec = poles_vec(find(abs(poles_vec) < 1));
P = length(poles_vec);

% Getting the number and length of chunks for the train set
ll_train = numel(chunk_train);   % number of chunks
L_Ch_train = zeros(ll_train,1);  %length of each chunk
for i=1:ll_train
    L_Ch_train(i) = length(chunk_train(i).stim);  
end

% Gamma for the train set
m_train = max (L_Ch_train);      % maximum chunk length
Gamma_train = zeros(m_train,P);
for row_index = 0:m_train-1
    Gamma_train(row_index+1,:) = poles_vec.^row_index;
end      
alpha_train = zeros(P,1);        % scaling factor
sGamma_train = zeros(m_train,P); % scaled Gamma
for pInd=1:P
    alpha_train(pInd) = (1-abs(poles_vec(pInd))^2)/(1-abs(poles_vec(pInd))^(2*mean(L_Ch_train))+2);
    sGamma_train(:,pInd) = Gamma_train(:,pInd)*alpha_train(pInd);
end

% For each chunk, creating structure arrays of X and sGamma corresponding to the input and output of that chunk
for l=1:ll_train
    X_train = tril(toeplitz(chunk_train(l).stim));   
    chunk_train(l).X = X_train;
    chunk_train(l).sGamma_train = sGamma_train(1:L_Ch_train(l),:);
end

% Getting the total output in double for the train set
output_train = [];
for kk= 1:numel(chunk_train)
    outputi_train = chunk_train(kk).out;
    output_train = [output_train; outputi_train];
end

% Finding the index of missing data in the output_train
idxi_train = find(output_train ~=1 & output_train ~=0);

% Solving system identification problem using L1 norm relaxation
y_hat_train = [];
cvx_begin
   variable c(P) complex                  % unknown input response coefficients of the system
   variable c_ic(P*ll_train) complex      % unknown natural response coefficients of the system
   variable t(P)
   minimize norm(t,1);
   subject to
   for l=1:ll_train
       yhat_train = real(chunk_train(l).X*(chunk_train(l).sGamma_train(1:end,:)*c) + chunk_train(l).sGamma_train*c_ic((l-1)*P+1:P*l,1));
       y_hat_train = [y_hat_train ; yhat_train];
   end
   for jj = 1:length(output_train)
       if jj ~= idxi_train
           (1+epsilon)*output_train(jj)+(1-output_train(jj))*(-1+epsilon)<=(y_hat_train(jj).*output_train(jj)-y_hat_train(jj).*(1-output_train(jj)))
       end
   end
   c(1:P_complex/2) == conj(c(P_complex/2+1:P_complex))
   for l=1:ll_train
       c_ic(1+P*(l-1):P*(l-1)+P_complex/2) == conj(c_ic(P*(l-1)+P_complex/2+1:P_complex+ P*(l-1)))
   end
   t2=repmat(t,ll_train,1);
   abs(c(1:P)) <= t;
   abs(c_ic) <= t2;
cvx_end
ind = find(abs(c)> prec);                 % finding significant coefficients (no less that the threshold "prec")
order = length(ind);

% Validation
% Getting the number and length of chunks for the test set
ll_test = numel(chunk_test);      % number of chunks
L_Ch_test = zeros(ll_test,1);     % length of each chunk
for i=1:ll_test
    L_Ch_test(i) = length(chunk_test(i).stim);  
end

% Gamma for the test set
m_test = max(L_Ch_test);          % maximum chunk length
Gamma_test = zeros(m_test,P);
for row_index = 0:m_test-1
    Gamma_test(row_index+1,:) = poles_vec.^row_index;
end      
alpha_test = zeros(P,1);          % scaling factor
sGamma_test = zeros(m_test,P);    % scaled Gamma
for pInd=1:P
    alpha_test(pInd) = (1-abs(poles_vec(pInd))^2)/(1-abs(poles_vec(pInd))^(2*mean(L_Ch_test))+2);
    sGamma_test(:,pInd) = Gamma_test(:,pInd)*alpha_test(pInd);
end

for l=1:ll_test
    X_test = tril(toeplitz(chunk_test(l).stim));   
    chunk_test(l).X_test = X_test;
    chunk_test(l).sGamma_test = sGamma_test(1:L_Ch_test(l),:);
end

% Getting the total number of samples
N = 0;
for i=1:ll_test
    N = N+length(chunk_test(i).stim);  
end
o = 0:N-1;

% Getting the total output in double for the test set
output_test = [];
for hh= 1:numel(chunk_test)
    outputi_test = chunk_test(hh).out;
    output_test = [output_test; outputi_test];
end

% Finding the index of missing data in the output_test
idxi_test = find(output_test ~=1 & output_test ~=0);

% System identification
y_hat_test = [];
i = 1;
j = L_Ch_test(1);
L_Ch_test(1+ll_test) = 0;

cvx_begin
   variable del(N)
   variable C_IC(P*(ll_test),1) complex        % unknown natural response coefficients of the system   
   minimize norm(del,1);
   subject to
   for l=1:ll_test
       yhat_test = real(chunk_test(l).X_test*(chunk_test(l).sGamma_test(1:end,:)*c) + chunk_test(l).sGamma_test*C_IC((l-1)*P+1:P*l,1)); 
       y_hat_test = [y_hat_test ; yhat_test];
       i = j+1;
       j = j+ L_Ch_test(l+1);
   end
   for ff = 1:length(output_test)
       if ff ~= idxi_test
           (1+epsilon)*output_test(ff)+(1-output_test(ff))*(-1+epsilon)<=((y_hat_test(ff)+del(ff)).*output_test(ff)-(y_hat_test(ff)+del(ff)).*(1-output_test(ff)))
       end
   end
   for l=1:ll_test
       C_IC(1+P*(l-1):P*(l-1)+P_complex/2) == conj(C_IC(P*(l-1)+P_complex/2+1:P_complex+ P*(l-1)))
   end 
cvx_end
error = norm(del,1);

% System 
transfer_function = 0;
for i = 1:order
transfer_function = transfer_function + c(ind(i))*tf([(1-abs(poles_vec(ind(i)))^2)/(1-abs(poles_vec(ind(i)))^(2*mean(L_Ch_train))+2)],[1 -poles_vec(ind(i))],Ts);
end
[num,den]=tfdata(transfer_function);
num = real(num{1});
den = real(den{1});
Lin_sys = tf(num,den,Ts);