function stuart_landau_delay_dynamics(n,para,sigma)
% This code simulates delay-coupled Stuart-Landau oscillators using MATLAB dde23
% n -- number of oscillators
% para -- which oscillator to make heterogeneous ('none', 'lambda', 'w', 'gamma', 'all')
% sigma -- level of heterogeneity measured by standard deviation
% For Fig.1 in the paper, n=18, para='w', sigma = 0.1.

% coupling parameters
beta = 0;
kmu = .3; k = .3; % coupling strength
tau = 1.8*pi; % coupling delay

% oscillator parameters
lambda0 = .1;
w0 = .12;
gamma0 = 0;
h = -.4;

% random seed
rng('shuffle')

% adjacency matrix for a directed ring of n nodes
adj = zeros(1,n);
adj(1,2) = 0; adj(1,end) = 1;
ring = adj;
for ii = 1:n-1
  adj = [adj;circshift(ring,[0,ii])];
end

% finding the sync solution
function y = f_omega(x)
y = w0 - gamma0*(lambda0 + kmu*(cos(beta-x*tau) - cos(beta)))...
	+ kmu*(sin(beta-x*tau) - sin(beta)) - x;
end

x0 = 1;
omega = fzero(@f_omega,x0) % frequency of the sync solution

phi = beta - omega*tau;
r2 = lambda0 + kmu*(cos(phi) - cos(beta)) % square amplitude of the sync solution

% adding heterogeneity to the specified parameter(s)
heterogeneity = normrnd(0,1,n,1);
heterogeneity = (heterogeneity - mean(heterogeneity))*sigma/std(heterogeneity);

lambda = lambda0*ones(n,1);
w = (w0+h)*ones(n,1);
gamma = (gamma0+h/r2)*ones(n,1);

switch para
	case 'lambda'
		lambda = lambda + heterogeneity
	case 'w'
		w = w + heterogeneity;
	case 'gamma'
		gamma = gamma + heterogeneity/r2;
	case 'all'
		lambda = lambda + heterogeneity;

		heterogeneity = normrnd(0,1,n,1);
		heterogeneity = (heterogeneity - mean(heterogeneity))*sigma/std(heterogeneity);
		w = w + heterogeneity;

		heterogeneity = normrnd(0,1,n,1);
		heterogeneity = (heterogeneity - mean(heterogeneity))*sigma/std(heterogeneity);
		gamma = gamma + heterogeneity/r2;
end

% simulating the system for T time units
T = 3000;
time = [0 T];
opts = ddeset('RelTol',1e-4,'AbsTol',1e-4);

% perturbation added to the sync state at t=0
perturb = 1e-1*rand(2*n,1);

% solving the DDE
sol = dde23(@rhs,tau,@rhshist,time,opts);

% extract phases of the oscillators
thetas = angle(sol.y(1:n,:)+i*sol.y(n+1:2*n,:));
% calculate Kuramoto order parameter
order1 = abs(sum(exp(i*thetas),1))/n;
% extract max amplitude
maxr = max(abs(sol.y(1:n,:)+i*sol.y(n+1:2*n,:)));
% calculate an order parameter that takes both phase and amplitude into account
order2 = abs(sum(sol.y(1:n,:)+i*sol.y(n+1:2*n,:),1))/n./maxr;
% time-averaged order parameters
sync1 = mean(order1(floor(end*.8):end))
sync2 = mean(order2(floor(end*.8):end))

% calculating the phase velocity $\hat{\Omega}_j$
for ii = 1:n
	phase_vel(ii,:) = (...
		thetas(ii,2:end) - thetas(ii,1:end-1)...
		+ 2*pi*(thetas(ii,2:end) - thetas(ii,1:end-1) < -pi)...
		- 2*pi*(thetas(ii,2:end) - thetas(ii,1:end-1) > pi)...
		)./(sol.x(2:end) - sol.x(1:end-1));
end


function dxdt = rhs(t,x,Z)  % set up DDEs
	dxdt = zeros(2*(n),1);
	xlag = Z(:,1);
	for ii = 1:n
		dxdt(ii) = lambda(ii)*x(ii) - w(ii)*x(n+ii) + ...
		(x(ii)^2 + x(n+ii)^2)*(gamma(ii)*x(n+ii) - x(ii)) + ...
		k*adj(ii,:)*(xlag(1:n) - x(ii));
		dxdt(n+ii) = lambda(ii)*x(n+ii) + w(ii)*x(ii) - ...
		(x(ii)^2 + x(n+ii)^2)*(gamma(ii)*x(ii) + x(n+ii)) + ...
		k*adj(ii,:)*(xlag(n+1:2*n) - x(n+ii));
	end
end

function S = rhshist(t)	% initial conditions
	histx = sqrt(r2)*cos(omega*t)*ones(n,1);
	histy = sqrt(r2)*sin(omega*t)*ones(n,1);
	S = [histx;histy] + perturb*heaviside(t);
end

% plotting the state trajectories, order parameters, and phase velocities of the oscillators
fonttype = 'Times';
fsize = 23;
txtattrib = {'FontName',fonttype,'FontSize',fsize,...
         'FontWeight','bold'};
txtattrib2 = {txtattrib{:},'Interpreter','Latex'};
set(0,'DefaultAxesFontSize',15)

figure

hAxis(1)=subplot(3,1,1);
pos = get( hAxis(1), 'Position');
pos(1)=.09;
pos(2)=.74;
pos(3)=.87;
pos(4)=.25;
set(hAxis(1), 'Position', pos);
hold on
for ii = 1:n
	plot(sol.x,sol.y(n+ii,:))
end
hold off
box on
ylabel('$y_j$', txtattrib2{:})
ylim([-.4,.4])
set(gca,'YTick',[-.3 0 .3])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gca,'xtick',[]);

hAxis(2)=subplot(3,1,2);
pos = get( hAxis(2), 'Position');
pos(1)=.09;
pos(2)=.46;
pos(3)=.87;
pos(4)=.25;
set(hAxis(2), 'Position', pos);
plot(sol.x,order1,'-b',sol.x,order2,'-r','LineWidth',2);
ylabel('$R_{1,2}$', txtattrib2{:})
ylim([0,1])
set(gca,'YTick',[0 .5 1])
%ylim([.8,1])
%set(gca,'YTick',[.8 .9 1])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gca,'xtick',[]);

hAxis(3)=subplot(3,1,3);
pos = get( hAxis(3), 'Position');
pos(1)=.09;
pos(2)=.18;
pos(3)=.87;
pos(4)=.25;
set(hAxis(3), 'Position', pos);
hold on
for ii = 1:n
	plot(sol.x(1:end-1),phase_vel(ii,:));
end
hold off
box on
xlabel('$t$', txtattrib2{:})
ylabel('$\Omega_j$', txtattrib2{:})
ylim([-2.5,2.5])
set(gca,'YTick',[-2 0 2])
%ylim([-.12,.12])
%set(gca,'YTick',[-.1 0 .1])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))

set(gcf, 'PaperPosition', [0 0 8 3]);
set(gcf, 'PaperSize', [8 3]);
saveas(gcf,strcat('trj_',para,'.pdf'));

end
