% PARAMETERS
nSim = 500;
x    = 10;
y    = 10;
roe  = 0.9;

% Set the seed.
randn('state', sum(100*clock));

samples = zeros(nSim,2);
sig     = 1 - roe^2;
n       = randn(2,nSim);

% Simulate the samples.
tic;
for i = 1:nSim,
  
  % Calculate the new x using the marginal (blocked move).
  x = n(1,i);
  
  % Calculate the new y using the conditional.
  y = roe*x + n(2,i) * sig;
  samples(i,:) = [x y];
end;

fprintf('Sampling took %0.3f seconds. \n', toc);

% Plot the values.
samples = samples';
plot(samples(1,:), samples(2,:), 'b.');
axis([-3 11 -3 11]);
xlabel('x');
ylabel('y');
