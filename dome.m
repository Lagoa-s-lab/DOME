function [Lin_sys] = dome (chunk,Ts,r,prec,epsilon)

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
   variable D
   variable c(P) complex                  % unknown input response coefficients of the system
   variable c_ic(P*ll) complex            % unknown natural response coefficients of the system
   variable t(P)
   minimize norm(t,1);
   subject to
   for l=1:ll
       yhat = real(chunk(l).X*(D+chunk(l).sGamma(1:end,:)*c) + chunk(l).sGamma*c_ic((l-1)*P+1:P*l,1));
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
transfer_function = D;
for i = 1:order
transfer_function = transfer_function + c(ind(i))*tf([(1-abs(poles_vec(ind(i)))^2)/(1-abs(poles_vec(ind(i)))^(2*mean(L_Ch))+2)],[1 -poles_vec(ind(i))],Ts);
end
[num,den]=tfdata(transfer_function);
num = real(num{1});
den = real(den{1});
Lin_sys = tf(num,den,Ts);
