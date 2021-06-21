function [Lin_sys] = dome (chunk,Ts,r,prec,epsilon)

% Creating poles consisting complex and real poles, with P being the number of poles
poles_real = linspace(-1,1,r);
poles_vec = [];
for rp = 1:r
    poles_vec = [poles_vec,poles_real(rp)+j*linspace(0.1,1,r)];
end
poles_vec = [poles_vec,conj(poles_vec)];
poles_vec = poles_vec(find(abs(poles_vec) < 1));
P_complex = length(poles_vec);
poles_vec = [poles_vec,linspace(-1,1,2*r)];
poles_vec = poles_vec(find(abs(poles_vec) < 1));
P = length(poles_vec);

% Getting the number and length of chunks
ll = numel(chunk);   % number of chunks
L_ch = zeros(ll,1);  %length of each chunk
for i=1:ll
    L_ch(i) = length(chunk(i).stim);  
end

% Gamma
m = max (L_ch);  % maximum chunk length
Gamma = zeros(m,P);
for row_index = 0:m-1
    Gamma(row_index+1,:) = poles_vec.^row_index;
end      
alpha = zeros(P,1);    % scaling factor
sGamma = zeros(m,P);   % scaled Gamma
for pInd=1:P
    alpha(pInd) = (1-abs(poles_vec(pInd))^2)/(1-abs(poles_vec(pInd))^(2*mean(L_ch))+2);
    sGamma(:,pInd) = Gamma(:,pInd)*alpha(pInd);
end

% Creating structure arrays of stimulus and outcome, and corresponding X and sGamma for each chunk
for l=1:ll
    X = tril(toeplitz(chunk(l).stim));   
    chunk(l).X = X;
    chunk(l).sGamma = sGamma(1:L_ch(l),:);
end

% Solving system identification problem using L1 norm relaxation 
cvx_begin
   variable D
   variable C(P*(ll+1),1) complex    % unknown coefficients of the system with the first P ones being the coefficients for the impulse response and the rest for the intrinsic response
   minimize norm(C,1);
   subject to
   for l=1:ll
       yhat = real(chunk(l).X*[D;chunk(l).sGamma(1:end-1,:)*C(1:P)] + chunk(l).sGamma*C(((l-1)*P+1:P*l)+P,1)); 
       (1+epsilon)*chunk(l).out+(1-chunk(l).out)*(-1+epsilon)<=(yhat.*chunk(l).out-yhat.*(1-chunk(l).out))
   end
   for l=1:ll+1
       C(1+P*(l-1):P*(l-1)+P_complex/2) == conj(C(P*(l-1)+P_complex/2+1:P_complex+ P*(l-1)))
   end 
cvx_end

ind = find(abs(C(1:P))> prec);  % finding significant coefficients (no less that the threshold "prec")
sp = length(ind);

% System 
transfer_function = D;
for i =1:sp
transfer_function = transfer_function + C(ind(i))*tf([(1-abs(poles_vec(ind(i)))^2)/(1-abs(poles_vec(ind(i)))^(2*mean(L_ch))+2)],[1 -poles_vec(ind(i))],Ts);
end

[num,den]=tfdata(transfer_function);
num = real(num{1});
den = real(den{1});
Lin_sys = tf(num,den,Ts);
