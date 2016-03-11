clear

n           = 10000;            % number of time points
dt          = .001;             % time setp
num_loops   = 100;              % repeat how many times
tau         = 2.3E-3;           % current time constant
alpha       = 1/.100;           % integration time constant

t           = (-1000:n-1)*dt;   % start time one second early 
f           = (0:n-1)/max(t);   % frequencies

% current response function
idx = t>=0 & t<.100;
sc  = t(idx).^.13.*exp((1-t(idx))/tau);

% output current
saved_I = zeros(n, num_loops);

for loop = 1:num_loops
    
    % input randomg noise
    Q  = randn(size(t));
    
    % convolve input with current response function
    Q = conv(Q, sc, 'full'); Q = Q(1:end-length(sc)+1);

    % dendritic current
    I  = zeros(size(t));
    I(1) = rand;
    
    % Loop over time (leaky integrator)
    for ii = 2:length(t)
        di = dt*(-alpha*I(ii-1)+Q(ii-1));
        I(ii) = I(ii-1) + di;
    end 
    
    % store this iteration
    saved_I(:, loop) = I(t>=0);
end

% Crop out pre trial data (t<0)
I = I(t>=0);
Q = Q(t>=0);
t = t(t>=0);

% Plot
figure(1); clf
subplot(2,2,1)
plot(t,I); title('I')

subplot(2,2,2)
plot(t,Q); title('Q')


subplot(2,2,3)
plot(f, abs(fft(Q))); title('Q')
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [.1 400])

subplot(2,2,4)
F_I = exp(mean(log(abs(fft(saved_I))),2));
plot(f, F_I);
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [.1 1000])
title('I')
hold on;
yl = get(gca, 'YLim');
plot(alpha/(2*pi)*[1 1], yl, 'r--')
plot(1/(tau*2*pi)*[1 1], yl, 'g--')
yl = get(gca, 'YLim');
text(alpha/(2*pi), yl(2)/2, 'Alpha', 'Color', 'r')
text(1/(tau*2*pi), yl(2)/2, 'Tau', 'Color', 'g')

idx =  f>100 & f<200;
x = log(f(idx));
y = log(F_I(idx));
mdl = fitlm(x, y);

pwr = 2*mdl.Coefficients.Estimate(2);
disp(pwr)


