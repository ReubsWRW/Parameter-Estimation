%% experimental data 
% ---------------------------------------------------------------
% THIS IS SAMPLE DATA, ACTUAL DATA IS TO BE INPUT HERE
x = -10:1:10;
noise = (0.1 - (-0.1)).*rand(1, length(x)) + (-0.1);
y = -5 * exp(-x.^2 / (2 * 3^2)) + noise;
% ---------------------------------------------------------------


%% Nonlinear Regression
% ---------------------------------------------------------------
% model ansatz (width)
fun = @(theta) -theta(1) * exp(-x.^2 / (2 * theta(2)^2)) - y;
options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-20);
init_guess = [1, 1]; % initial guess of parameters
theta = lsqnonlin(fun, init_guess, [], [], options); % nonlinear least-square regression

% printing computed parameters
disp('GAUSSIAN PARAMETERS: ')
h = theta(1); % Gaussian height
disp(['h: ', num2str(theta(1))]) 
sigma = theta(2); % Gaussian width
disp(['sigma: ', num2str(theta(2))]) 

% computing model plot from computed parameters
dx = 0.05; % numerical step size
model_x = x(1):dx:x(end);
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
x_r0 = [2.25, 2.5]; % initial guess of parameters (HAS TO BE TUNED!!)
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
% single-Gaussian channel
subplot(1, 2, 1) 
plot(model_x, model, 'b')

hold on; % plotting actual data points
scatter(x, y, 'ro')
hold on; % plotting largest particles
circle(0, -r_sg , r_sg, 'g-') % single-Gaussian largest particle

% axis labels and limits
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
xlim([x(1) x(end)])
ylim([x(1) x(end)])


subplot(1, 2, 2)
% dual-Gaussian channel
plot(model_x, model, 'b')
hold on;
plot(model_x, -model, 'b')

hold on; % plotting actual data points
scatter(x, y, 'ro')
hold on;
circle(0, 0, r, 'g-'); % dual-Gaussian largest particle

% axis labels and limits
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
xlim([x(1) x(end)])
ylim([x(1) x(end)])
% ---------------------------------------------------------------


%% additional function definitions
% ---------------------------------------------------------------
function circle(x,y,r,color)
    % Description: 
    %   function that plots a circle given input parameters.
    %
    % Arguments: 
    %   - x, y: coordinates to define center of the circle.
    %   - r: radius of the circle. 
    %   - color: to specify plotting color of circle
    
    ang=0:0.01:2*pi; 
    xp = r * cos(ang);
    yp = r * sin(ang);
    plot(x+xp,y+yp,color);
end
% ---------------------------------------------------------------
    

