% Linearization of the Cortical Multistable Neuron Model with Phase Portrait and Time Plots

clear; clc; close all;
disp('--------------------------------------------------------------')
warning('off')

% Time -----
t0 = 0; tmax = 10; dt = 1e-3; % s
T = (t0:dt:tmax)'; nT = length(T);

% Symbolic variables -----
syms x1_ x2_
vars_ = [x1_ x2_];

% Model -----
Model_Name = 'Cortical Multistable Neuron';
par = []; par_ = [];

% Initial conditions
x10 = 0.1; x20 = 0.01;
xmin = 0; xmax = 60;

% Model equations
f1 = (-x1_ + 50/(1 + exp(-(x2_*x1_ - 5)/2)))/10;
f2 = 0.02*x1_*(1 - x2_) - 0.01*x2_;
dxdt_ = [f1; f2];

% Display model
disp(['*** ' Model_Name ' ***'])
disp(['Function: dx1dt = ' char(simplify(dxdt_(1)))])
disp(['Function: dx2dt = ' char(simplify(dxdt_(2)))]);

% Solve fixed points
FP_ = solve(dxdt_ == 0, vars_);
Jacobian_ = jacobian(dxdt_, vars_);

% Initialize solution arrays
x1 = zeros(1,nT); x2 = zeros(1,nT);
x1(1) = x10; x2(1) = x20;

% Evaluate Jacobian, eigenvalues/vectors
nFP = length(FP_.x1_);
Stability = zeros(1,nFP);
Nature = cell(1,nFP);
for kFP = 1:nFP
    J_FP_ = subs(Jacobian_, vars_, [FP_.x1_(kFP), FP_.x2_(kFP)]);
    [eigvect, eigenval] = eig(J_FP_);
    eigs = diag(eigenval);
    if isreal(eigs)
        if eigs(1)>0 && eigs(2)>0; Nature{kFP} = 'unstable node'; Stability(kFP) = 1; end
        if sign(eigs(1)*eigs(2))<=0; Nature{kFP} = 'saddle node'; Stability(kFP) = 2; end
        if eigs(1)<=0 && eigs(2)<=0; Nature{kFP} = 'stable node'; Stability(kFP) = 3; end
    else
        if real(eigs(1))>0; Nature{kFP} = 'unstable focus'; Stability(kFP) = 4; end
        if real(eigs(1))==0; Nature{kFP} = 'center'; Stability(kFP) = 5; end
        if real(eigs(1))<0; Nature{kFP} = 'stable focus'; Stability(kFP) = 6; end
    end
end

% Numerical simulation
f1_num = matlabFunction(f1, 'Vars', {x1_, x2_});
f2_num = matlabFunction(f2, 'Vars', {x1_, x2_});
for kt = 1:nT-1
    x1(kt+1) = x1(kt) + dt * f1_num(x1(kt), x2(kt));
    x2(kt+1) = x2(kt) + dt * f2_num(x1(kt), x2(kt));
end

% Graphics settings
FontSize = 20; FP_FontSize = 18; LineWidth = 3; FP_Size = 20; IC_Size = 10;
FP_FaceColor = ['w' 'w' 'r' 'w' 'w' 'r']; FP_Symbol = [' ' 'x' ' ' ' ' '.' ' '];
Color_Map = turbo(256);
x_All = linspace(xmin, xmax, 100);
y_All = linspace(0, 1, 100);
[x1_Grad,x2_Grad] = meshgrid(x_All, y_All);
dx1dt_Grad = f1_num(x1_Grad,x2_Grad);
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
title(['Phase Portrait: ' Model_Name])
image('XData',x_All,'YData',y_All,'CData',Field_Color,'CDataMapping','direct')
set(gca,'YDir','normal'); hold on; box on;
fimplicit(f1 == 0, [xmin xmax 0 1],'c','LineWidth',LineWidth);
fimplicit(f2 == 0, [xmin xmax 0 1],'g','LineWidth',LineWidth);
plot(x1, x2, 'b', 'LineWidth', LineWidth);
plot(x1(1), x2(1), 'bo', 'MarkerSize', IC_Size, 'MarkerFaceColor', 'b');

% Fixed points
FP.x1 = double(FP_.x1_); FP.x2 = double(FP_.x2_);
for kFP = 1:nFP
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
