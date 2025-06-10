% Linearization of the Cortical Multistable Neuron Model with Phase Portrait and Time Plots

clear; clc; close all;
disp('--------------------------------------------------------------')
warning('off')

% Time -----
t0 = 0; tmax = 10; dt = 1e-3; % s
T = (t0:dt:tmax)'; nT = length(T);

% Plasticity mode ('anti' for Egorov-style, 'homeo' for depressive)
plasticity_type = 'homeo';

% Symbolic variables -----
syms x1_ x2_ input_
vars_ = [x1_ x2_];

% Model -----
Model_Name = 'Cortical Multistable Neuron';
par = []; par_ = [];

% Initial conditions
x10 = 0.1; x20 = 0.01;
xmin = 0; xmax = 60;

% Plasticity rule
if strcmp(plasticity_type, 'anti')
    f2 = 0.02*x1_*(1 - x2_) - 0.01*x2_;
elseif strcmp(plasticity_type, 'homeo')
    f2 = -0.02*x1_*(1 - x2_) - 0.01*x2_;
else
    error('Unknown plasticity type');
end

% Model equations with external input term
f1 = (-x1_ + 50/(1 + exp(-(x2_*x1_ - 5)/2)) + input_) / 10;
dxdt_ = [f1; f2];

% Display model
disp(['*** ' Model_Name ' ***'])
disp(['Function: dx1dt = ' char(simplify(f1))])
disp(['Function: dx2dt = ' char(simplify(f2))]);

% Solve fixed points (without input)
FP_ = solve(subs(dxdt_, input_, 0) == 0, vars_);
nFP = length(FP_.x1_);
Jacobian_ = jacobian(subs(dxdt_, input_, 0), vars_);

% Initialize solution arrays
x1 = zeros(1,nT); x2 = zeros(1,nT);
x1(1) = x10; x2(1) = x20;

% Phasic input protocol
pulse_duration = 0.25; % s
pulse_interval = 1.0;  % s
pulse_times = 0.5:pulse_interval:(tmax - pulse_duration);

% Numerical simulation
f1_input = matlabFunction(f1, 'Vars', {x1_, x2_, input_});
f2_num = matlabFunction(f2, 'Vars', {x1_, x2_});
for k_t = 1:nT-1
    t_curr = T(k_t);
    input = 0;
    if any(t_curr >= pulse_times & t_curr < pulse_times + pulse_duration)
        input = 5; % pulse amplitude
    end
    dx1 = f1_input(x1(k_t), x2(k_t), input);
    dx2 = f2_num(x1(k_t), x2(k_t));
    x1(k_t+1) = x1(k_t) + dt * dx1;
    x2(k_t+1) = x2(k_t) + dt * dx2;
end

% Graphics settings
FontSize = 20; FP_FontSize = 18; LineWidth = 3; FP_Size = 20; IC_Size = 10;
FP_FaceColor = ['w' 'w' 'r' 'w' 'w' 'r']; FP_Symbol = [' ' 'x' ' ' ' ' '.' ' '];
Color_Map = turbo(256);
x_All = linspace(xmin, xmax, 100);
y_All = linspace(0, 1, 100);
[x1_Grad,x2_Grad] = meshgrid(x_All, y_All);
dx1dt_Grad = f1_input(x1_Grad, x2_Grad, 0);
dx2dt_Grad = f2_num(x1_Grad,x2_Grad);
dpdt_theta = angle(dx1dt_Grad + 1i*dx2dt_Grad);
dpdt_theta_n = ceil((255)*(dpdt_theta + pi)/(2*pi)) + 1;
Field_Color = zeros([size(dpdt_theta_n), 3]);
for k = 1:3
    temp = Color_Map(:,k);
    Field_Color(:,:,k) = reshape(temp(dpdt_theta_n), size(dpdt_theta_n));
end

% Phase Portrait with color-coded direction field
figure(1); clf; set(gcf,'Color','white'); set(gca,'FontSize',FontSize)
title(['Phase Portrait: ' Model_Name ' (' plasticity_type ')'])
image('XData',x_All,'YData',y_All,'CData',Field_Color,'CDataMapping','direct')
set(gca,'YDir','normal'); hold on; box on;
fimplicit(subs(f1, input_, 0) == 0, [xmin xmax 0 1],'c','LineWidth',LineWidth);
fimplicit(f2 == 0, [xmin xmax 0 1],'g','LineWidth',LineWidth);
plot(x1, x2, 'b', 'LineWidth', LineWidth);
plot(x1(1), x2(1), 'bo', 'MarkerSize', IC_Size, 'MarkerFaceColor', 'b');

% Fixed points
FP.x1 = double(FP_.x1_); FP.x2 = double(FP_.x2_);
for kFP = 1:nFP
    J_FP_ = subs(Jacobian_, vars_, [FP.x1(kFP), FP.x2(kFP)]);
    [~, D] = eig(double(J_FP_)); eigs = diag(D);
    if isreal(eigs)
        if eigs(1)>0 && eigs(2)>0; Stability(kFP) = 1; end
        if sign(eigs(1)*eigs(2))<=0; Stability(kFP) = 2; end
        if eigs(1)<=0 && eigs(2)<=0; Stability(kFP) = 3; end
    else
        if real(eigs(1))>0; Stability(kFP) = 4; end
        if real(eigs(1))==0; Stability(kFP) = 5; end
        if real(eigs(1))<0; Stability(kFP) = 6; end
    end
    plot(FP.x1(kFP), FP.x2(kFP), 'ro', 'MarkerSize', FP_Size, 'LineWidth', LineWidth, 'MarkerFaceColor', FP_FaceColor(Stability(kFP)));
    plot(FP.x1(kFP), FP.x2(kFP), ['r' FP_Symbol(Stability(kFP))], 'LineWidth', LineWidth, 'MarkerSize', FP_Size);
end
xlabel('r (Hz)'); ylabel('w'); axis([xmin xmax 0 1]); axis square; grid on;

% Time series with fixed point lines
figure(2); clf; set(gcf,'Color','white');
subplot(2,1,1); hold on; box on; set(gca,'FontSize',FontSize)
title('Time Evolution: r');
plot(T, x1, 'r', 'LineWidth', 2);
for kFP = 1:nFP
    plot([0 tmax], [FP.x1(kFP) FP.x1(kFP)], '--k', 'LineWidth', 1.5)
end
ylabel('r (Hz)'); grid on;

subplot(2,1,2); hold on; box on; set(gca,'FontSize',FontSize)
title('Time Evolution: w');
plot(T, x2, 'b', 'LineWidth', 2);
for kFP = 1:nFP
    plot([0 tmax], [FP.x2(kFP) FP.x2(kFP)], '--k', 'LineWidth', 1.5)
end
xlabel('time (s)'); ylabel('w'); grid on;

disp('--------------------------------------------------------------')
warning('on')
