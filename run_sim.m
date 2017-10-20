% Tahmid Mehdi
% runs a simulation of the Master-Slave Latch for a simplified Drosophila Gap Gene Network
% April 2016

% initialize parameters
global kA2; global kKA; global lambdap; global kB2; global kJB; global kJ2; global kK2;
global v; global Khomodimers; global Kheterodimers; global lambda;
kA2 = [0.2 2]; kKA=[0.2 2]; lambdap=0.138; kB2=[0.2 2]; kJB=[0.2 2];
kJ2=[0.2 2]; kK2=[0.2 2]; v=[11.6 11.6 11.6 11.6]; Khomodimers=[8.3 7.3 8.1 7.3]; Kheterodimers=[7.65 7.65];
lambda=0.18;
% specification of dynamics
ODEFUN = @mslatch;

% simulation
[t,G]=ode45(ODEFUN, [0,130], [5,50,0,0]);

% Time-Concentration
fig1=figure(1)
plot(t,G(:,1),'m', t,G(:,2),'r', t,G(:,3),'c--', t,G(:,4),'g--', 'LineWidth', 2)
axis([0 130 0 50])
xlabel('Time (min)')
ylabel('Concentration (nM)')
legend('Hunchback', 'Knirps', 'Kruppel', 'Giant')

% Phase Portraits
fig1=figure(2)
plot(G(:,1), G(:,2), 'LineWidth', 2)
axis([0 50 0 50])
xlabel('[Hunchback] (nM)')
ylabel('[Knirps] (nM)')

fig1=figure(3)
plot(G(:,3), G(:,4), 'LineWidth', 2)
axis([0 50 0 50])
xlabel('[Kruppel] (nM)')
ylabel('[Giant] (nM)')
