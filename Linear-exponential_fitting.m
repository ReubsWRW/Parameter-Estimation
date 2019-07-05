%% experimental data 
% ---------------------------------------------------------------
width = [186 149 135 114 104 100 106;
        204	170	151	127	134	120	120;
        206	172	160	131	139	140	141;
        216	186	180	160	158	146	145;
        225	210	206	184	178	171	177]; % widt matrix
    
depth = [250 150 91	76 67 68 65;
        342	190	143	79	84	109	108;
        404	256	152	127	136	130	147;
        453	294	211	138	136	142	146;
        576	485	335	227	173	183	179]; % height matrix
    
p_data = [40; 50; 60; 70; 80]; % vector of powers used
v_data = [20; 30; 40; 50; 60; 70; 80]; % vector of speeds used
% ---------------------------------------------------------------


%% vetorizing data 
% -------------------------------------------------------------
v_data_new = [];
p_data_new = [];
for i = 1:length(v_data)
    for j = 1:length(p_data)
        v_data_new = [v_data_new; v_data(i)];
    end
end
for j = 1:length(v_data)
    p_data_new = [p_data_new; p_data];
end
x = [p_data_new, v_data_new]; % vectorized power/speed data

% vectorizing width values
y_w = []; % declaring vectorized width values
shape = size(width);
for col = 1:shape(2)
   y_w = [y_w; width(:, col)];
end

% vectorizing depth values
y_h = []; % declaring vectorized depth values
shape = size(depth);
for col = 1:shape(2)
   y_h = [y_h; depth(:, col)];
end
% -------------------------------------------------------------

%% Nonlinear Regression 
% ----------------------------------------------------------------
% model ansatz (width)
fun1 = @(t) (t(1)*x(:,1) + t(2) * ones(length(x(:,1)),1)) .* ...
            (exp(-t(3)*x(:,2)) + t(4)*ones(length(x(:,2)),1)) ...
            - y_w;
options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-20);
theta = lsqnonlin(fun1, [1, 1, 1, 1], [], [], options); % nonlinear least-square regression

% printing computed parameters
disp('WIDTH PARAMETERS: ')
disp(['alpha1: ', num2str(theta(1))])
disp(['alpha2: ', num2str(theta(2))])
disp(['beta: ', num2str(theta(3))])
disp(['gamma: ', num2str(theta(4))])

% paramters from nonlinear regression of width
alpha1 = theta(1);
alpha2 = theta(2);
beta = theta(3);
gamma = theta(4);

% computing model plot from computed parameters
model_w = zeros(length(p_data), length(v_data));
for row = 1:length(p_data)
   for col = 1:length(v_data)
       
       % populating width model matrix
       model_w(row, col) = (alpha1*p_data(row) + alpha2) * ...
                           (exp(- beta*v_data(col)) + gamma);
   
   end
end


% model ansatz (depth)
fun2 = @(t) (t(1)*x(:,1) + t(2)*ones(length(x(:,1)),1)) .* ...
            (exp(-t(3)*x(:,2)) + t(4)*ones(length(x(:,2)),1)) ...
            - y_h;
options = optimoptions('lsqnonlin', 'FunctionTolerance', 1e-20);
theta = lsqnonlin(fun2, [2.0976, -0.43455, 0.032866, 1.4281], [], [], options); % nonlinear least-square regression

% printing computed parameters
disp('HEIGHT PARAMETERS: ')
disp(['alpha1: ', num2str(theta(1))])
disp(['alpha2: ', num2str(theta(2))])
disp(['beta: ', num2str(theta(3))])
disp(['gamma: ', num2str(theta(4))])

% paramters from nonlinear regression of width
alpha1 = theta(1);
alpha2 = theta(2);
beta = theta(3);
gamma = theta(4);

% computing model plot from computed parameters
model_h = zeros(length(p_data), length(v_data));
for row = 1:length(p_data)
   for col = 1:length(v_data)
       
       % populating depth model matrix
       model_h(row, col) = (alpha1*p_data(row) + alpha2) * ...
                           (exp(- beta*v_data(col)) + gamma);
       
   end
end
% ----------------------------------------------------------------

%% plotting 
% ------------------------------------------------------------------
figure;
[X,Y] = meshgrid(20:10:80, 40:10:80); 

% SUBPLOT 1
% -----------------------------------------------
subplot(2,2,1)
surf(X, Y, width) % surface plot of width data

% axis labels
xlabel('speed')
ylabel('power')
title('Width Data')
colorbar
% -----------------------------------------------

% SUBPLOT 2
% -----------------------------------------------
subplot(2, 2, 2)
scatter3 (x(:,2), x(:,1), y_w, 'o', 'MarkerFaceColor', 'r'); % scatter plot od data
hold on
surf(X, Y, model_w) % surface plot of width model

% axis labels
xlabel('speed')
ylabel('power')
title('Width Model')
colorbar
% -----------------------------------------------

% SUBPLOT 3
% -----------------------------------------------
subplot(2,2,3)
surf(X, Y, depth) % surface plot of height data

% axis labels
xlabel('speed')
ylabel('power')
title('Height Data')
colorbar
% -----------------------------------------------

% SUBPLOT 4
% -----------------------------------------------
subplot(2,2,4)
scatter3 (x(:,2), x(:,1), y_h, 'o', 'MarkerFaceColor', 'r'); % scatter plot od data
hold on 
surf(X, Y, model_h) % surface plot of height model

% axis labels
xlabel('speed')
ylabel('power')
title('Height Model')
colorbar
% -----------------------------------------------
% ------------------------------------------------------------------











