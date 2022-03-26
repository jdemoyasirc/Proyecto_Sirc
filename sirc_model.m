clc;
clear;
close all;

% Model parameters
alfa = 1/14; %Average time spent by subjects in the compartment I 

delta = 1/180; %Average time spent by subjects in the compartment R
gamma = 1/14; %Average time spent by subjects in the compartment C

mu = 0.06; %Mortality rate in every compartment, equal to rate of newborn in population
r = 0.026; % rate of infection
Ro = 2.5; %Average reproduction number

beta= Ro *(mu + alfa); %Contact rate
sigma = r * (1/delta + 1/gamma); %Average reinfection probability of a cross-inmune subject

N = 1; % Total population N = S + E + I + R
I0 = 0.02; % initial number of infected
T = 30; % period of 30 days
dt = 1; % time interval of 6 hours (1/4 of a day)

fprintf('Value of alfa is %.2f \n',alfa)
fprintf('Value of delta is %.2f \n',delta)
fprintf('Value of gamma is %.2f \n',gamma)
fprintf('Value of mu is %.2f \n',mu)
fprintf('Value of sigma is %.2f \n',sigma)
fprintf('Value of beta is %.2f \n',beta)

fprintf('Value of parameter R0 calculated is %.2f \n',N*beta/gamma)
fprintf('Value of parameter R0 is %.2f \n',Ro)
fprintf('Value of parameter I0 is %.2f \n',I0)

% Calculate the model
[S,I,R,C] = sir_model(alfa, delta, gamma, beta, mu, sigma,N,I0,T,dt);
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,S,'b',tt,I,'r',tt,R,'g',tt,C,'y','LineWidth',1); %grid on;
xlabel('Days'); ylabel('Number of individuals');
legend('S','I','R','C');

figure
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,S,'b','LineWidth',1); %grid on;
xlabel('Days'); ylabel('Number of Susceptible');
legend('S');

figure
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,I,'r','LineWidth',1); %grid on;
xlabel('Days'); ylabel('Number of Infected');
legend('I');

figure
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,R,'g','LineWidth',1); %grid on;
xlabel('Days'); ylabel('Number of Recovered');
legend('R');

figure
% Plots that display the epidemic outbreak
tt = 0:dt:T-dt;
% Curve
plot(tt,C,'y','LineWidth',1); grid on;
xlabel('Days'); ylabel('Number of Crossed Inmunity');
legend('C');

function [S,I,R,C] = sir_model(alfa, delta, gamma, beta, mu, sigma,N,I0,T,dt)
    % if delta = 0 we assume a model without immunity loss
    S = zeros(1,T/dt);
    S(1) = N;
    I = zeros(1,T/dt);
    I(1) = I0;
    R = zeros(1,T/dt);
    R(1)=0.0;
    C = zeros(1,T/dt);
    C(1)=0.0;
    fprintf('tt | S | dS | I | dI | R | dR | C | dC \n')    
    for tt = 1:(T/dt)-1
        % Equations of the model
        dS = (mu * (1 - S(tt)) - beta * S(tt) * I(tt) + gamma * C(tt)) * dt;

        dI = (beta * S(tt) * I(tt) + sigma * beta * C(tt) * I(tt) - (mu + alfa) *  I(tt)) * dt;

        dR = ((1 - sigma) * beta * C(tt) * I(tt)  + alfa * I(tt) - (mu + delta) * R(tt)) * dt;

        dC = (delta * R(tt) - beta * C(tt) * I(tt) - (mu + gamma) * C(tt)) * dt;

        S(tt+1) = S(tt) + dS;
        I(tt+1) = I(tt) + dI;
        R(tt+1) = R(tt) + dR;
        C(tt+1) = C(tt) + dC;
        fprintf('%d | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f | %.2f\n',tt,S(tt),dS,I(tt),dI,R(tt),dR,C(tt),dC)
    end
end