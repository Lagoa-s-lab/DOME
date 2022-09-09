function [Lin_sys] = dome (chunk,Ts,r,prec,epsilon)

% This function takes 5 inputs, chunk,Ts,r,prec,epsilon, and returns Lin_sys.
% Lin_sys = dome (chunk,Ts,r,prec,epsilon) returns the transfer function representation of the linear system identified by the algorithm.
% The structure of inputs to the "dome" function is as follows: 
% 1. "chunk" is a structured array whose size is the number of chunks and where each of the entries has two elements in its structure,
% "chunk(i).stim" and "chunk(i).out". These are two arrays of the same size with the chunk(i).stim containing the values of the
% stimulus for the i-th chunk and chunk(i).out containing the corresponding outcomes.
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

% Getting the number and length of chunks
ll = numel(chunk);   % number of chunks
L_Ch = zeros(ll,1);  %length of each chunk
for i=1:ll
    L_Ch(i) = length(chunk(i).stim);  
end

% Gamma
m = max (L_Ch);  % maximum chunk length
Gamma = zeros(m,P);
for row_index = 0:m-1
    Gamma(row_index+1,:) = poles_vec.^row_index;
end      
alpha = zeros(P,1);    % scaling factor
sGamma = zeros(m,P);   % scaled Gamma
for pInd=1:P
    alpha(pInd) = (1-abs(poles_vec(pInd))^2)/(1-abs(poles_vec(pInd))^(2*mean(L_Ch))+2);
    sGamma(:,pInd) = Gamma(:,pInd)*alpha(pInd);
end

% For each chunk, creating structure arrays of X and sGamma corresponding to the input and output of that chunk
for l=1:ll
    X = tril(toeplitz(chunk(l).stim));   
    chunk(l).X = X;
    chunk(l).sGamma = sGamma(1:L_Ch(l),:);
end

% Getting the total output in double
output = [];
for kk= 1:numel(chunk)
    outputi = chunk(kk).out;
    output = [output; outputi];
end

% Finding the index of missing data in the output
idxi = find(output ~=1 & output ~=0);

% Solving system identification problem using L1 norm relaxation 
y_hat = [];
cvx_begin
   variable c(P) complex                  % unknown input response coefficients of the system
   variable c_ic(P*ll) complex            % unknown natural response coefficients of the system
   variable t(P)
   minimize norm(t,1);
   subject to
   for l=1:ll
       yhat = real(chunk(l).X*(chunk(l).sGamma(1:end,:)*c) + chunk(l).sGamma*c_ic((l-1)*P+1:P*l,1));
       y_hat = [y_hat ; yhat];            
   end
   for jj = 1:length(output)
       if jj ~= idxi
           (1+epsilon)*output(jj)+(1-output(jj))*(-1+epsilon)<=(y_hat(jj).*output(jj)-y_hat(jj).*(1-output(jj)))
       end
   end
   c(1:P_complex/2) == conj(c(P_complex/2+1:P_complex))
   for l=1:ll
       c_ic(1+P*(l-1):P*(l-1)+P_complex/2) == conj(c_ic(P*(l-1)+P_complex/2+1:P_complex+ P*(l-1)))
   end   
   t2=repmat(t,ll,1);
   abs(c(1:P)) <= t;
   abs(c_ic) <= t2;
cvx_end
ind = find(abs(c)> prec);                 % finding significant coefficients (no less that the threshold "prec")
order = length(ind);

% System 
transfer_function = 0;
for i = 1:order
transfer_function = transfer_function + c(ind(i))*tf([(1-abs(poles_vec(ind(i)))^2)/(1-abs(poles_vec(ind(i)))^(2*mean(L_Ch))+2)],[1 -poles_vec(ind(i))],Ts);
end
[num,den]=tfdata(transfer_function);
num = real(num{1});
den = real(den{1});
Lin_sys = tf(num,den,Ts);
