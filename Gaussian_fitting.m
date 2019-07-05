%% experimental data 
% ---------------------------------------------------------------
% THIS IS FAKE DATA, INPUT YOUR ACTUAL DATA HERE
x = -10:1:10;
% y = [];
noise = (0.1 - (-0.1)).*rand(1, length(x)) + -0.1;
y = -5 * exp(-x.^2 / (2 * 3^2)) + noise;
% ---------------------------------------------------------------


%% Nonlinear Regression
% ---------------------------------------------------------------
% model ansatz (width)
fun = @(theta) -theta(1) * exp(-x.^2 / (2 * theta(2)^2)) - y;
options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-20);
theta = lsqnonlin(fun, [1, 1], [], [], options); % nonlinear least-square regression

% printing computed parameters
disp('GAUSSIAN PARAMETERS: ')
h = theta(1); % Gaussian height
disp(['h: ', num2str(theta(1))]) 
sigma = theta(2); % Gaussian width
disp(['sigma: ', num2str(theta(2))]) 

% computing model plot from computed parameters
model_x = x(1):0.05:x(end);
model = [];
for pts = 1:length(model_x)
    model = [model, -h * exp(-model_x(pts).^2 / (2 * sigma^2))];
end
% ---------------------------------------------------------------

%% largest particle radius (dual-Gaussian)
% ---------------------------------------------------------------
if (h >= sigma)
   r = sigma * sqrt(1 + 2*log(h/sigma));
else
   r = h;
end

disp (['Largest particle radius (dual-Gaussian): ', num2str(r)]);
% ---------------------------------------------------------------


%% largest particle radius (single-Gaussian)
% ---------------------------------------------------------------
F = @(x_r) x_r(2) + sqrt(x_r(2)^2 - x_r(1)^2) - h*exp(-x_r(1)^2/(2*sigma^2));
options = optimset('Display','off', 'TolX', 1e-20);
x_r0 = [2.25, 2.5]; % initial conditions have to be tuned!!
x_r = fminsearch(F, x_r0, options);
if (h < sigma)
    r_sg = h/2;    
else
    r_sg = x_r(2);
end
disp(['Largest particle radius (single-Gaussian): ', num2str(r_sg)])
% ---------------------------------------------------------------


%% Plotting
% ---------------------------------------------------------------
figure;
plot(model_x, model)
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')

hold on

scatter(x, y, 'ro')

hold on 

circle(0, 0, r); % dual-Gaussian largest particle

hold on

circle(0, -r_sg , r_sg) % single-Gaussian largest particle
% ---------------------------------------------------------------


function circle(x,y,r)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp = r * cos(ang);
    yp = r * sin(ang);
    plot(x+xp,y+yp, 'g');
end
    

