%% Plot sample activation functions
% last modified: June 10, 2020

x = -10:.01:10;
figure
% plot sign
subplot(2,2,1)
plot(x,sign(x),'-')
ax = gca;
ax.FontSize = 15;
yticks([-1 0 1])
xticks([-5 0 5])
title("heaviside", 'fontname', 'calibri')

% plot hyperbolic tan
subplot(2,2,2)
plot(x,tanh(x))
ax = gca;
ax.FontSize = 15;
yticks([-1 0 1])
xticks([-5 0 5])
title("hyperbolic tangent", 'fontname', 'calibri')

% plot sigmoid
subplot(2,2,3)
plot(x,sigmoid(x))
ax = gca;
ax.FontSize = 15;
yticks([-1 0 1])
xticks([-5 0 5])
title("sigmoid", 'fontname', 'calibri')

% plot reLU
impulse = x==0;
unitstep = x>=0;
ramp = x.*unitstep;
quad = x.^2.*unitstep;

subplot(2,2,4)
plot(x,ramp)
ax = gca;
ax.FontSize = 15;
yticks([-1 0 1])
xticks([-5 0 5])
title("reLU", 'fontname', 'calibri')