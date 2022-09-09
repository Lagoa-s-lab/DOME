%% Plots

% Impulse Response
figure (1)
[y,t] = impulse(Lin_sys);
plot(t,y,'color','[0.9290, 0.6940, 0.1250]','LineWidth', 2)
title('System Impulse Response');xlabel('Time (minute)'); ylabel('Impulse Response');

% Pulse Response
figure (2)
stim = [ones(10,1); zeros(length(y)-10, 1)];  % considering 10 pulses    
h = tril(toeplitz(y));
Pulsei= h*stim;
plot (t, Pulsei,'color','[0.9290, 0.6940, 0.1250]','LineWidth', 2);
title('System Pulse Response');
xlabel('Time (minute)'); ylabel('Pulse Response');