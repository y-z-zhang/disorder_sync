function sl_delay_single(n,para)
% This code simulates the delay-coupled SL oscillators using MATLAB dde23

fonttype = 'Times';
fsize = 23;
txtattrib = {'FontName',fonttype,'FontSize',fsize,...
         'FontWeight','bold'};
txtattrib2 = {txtattrib{:},'Interpreter','Latex'};

% parameters
beta = 0;
kmu = .3; k = .3;
tau = 1.8*pi;

lambda0 = .1;
w0 = .12;
%w0 = 1;
gamma0 = 0;
h = -.4;
%h = .9;
hetero = .2;

rng('shuffle')

% adjaency matrix
adj = zeros(1,n);
adj(1,2) = 0; adj(1,end) = 1;
ring = adj;
for ii = 1:n-1
  adj = [adj;circshift(ring,[0,ii])];
end

function y = f_omega(x)
%% solve omega as x
y = w0 - gamma0*(lambda0 + kmu*(cos(beta-x*tau) - cos(beta)))...
	+ kmu*(sin(beta-x*tau) - sin(beta)) - x;
end

x0 = 1;
omega = fzero(@f_omega,x0)

phi = beta - omega*tau;
r2 = lambda0 + kmu*(cos(phi) - cos(beta))

R = [cos(phi) -sin(phi); sin(phi) cos(phi)];


sample = normrnd(0,1,n,1);
sample = (sample - mean(sample))*hetero/std(sample)
lambda = lambda0*ones(n,1);
w = (w0+h)*ones(n,1);
gamma = (gamma0+h/r2)*ones(n,1);
switch para
	case 'lambda'
		lambda = lambda + sample
	case 'w'
		w = w + sample;
	case 'gamma'
		gamma = gamma + sample/r2;
	case 'all'
		lambda = lambda + sample;
		sample = normrnd(0,1,n,1);
		sample = (sample - mean(sample))*hetero/std(sample);
		w = w + sample;
		sample = normrnd(0,1,n,1);
		sample = (sample - mean(sample))*hetero/std(sample);
		gamma = gamma + sample/r2;
end

T = 300;
time = [0 T];
%opts = ddeset('RelTol',1e-3,'AbsTol',1e-6);
opts = ddeset('RelTol',1e-4,'AbsTol',1e-4);

amplitude = 1e-1*rand(2*n,1);

sol = dde23(@rhs,tau,@rhshist,time,opts);

thetas = angle(sol.y(1:n,:)+i*sol.y(n+1:2*n,:));
order1 = abs(sum(exp(i*thetas),1))/n;
maxr = max(abs(sol.y(1:n,:)+i*sol.y(n+1:2*n,:)));
order2 = abs(sum(sol.y(1:n,:)+i*sol.y(n+1:2*n,:),1))/n./maxr;
sync1 = mean(order1(floor(end*.8):end))
sync2 = mean(order2(floor(end*.8):end))

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
  		%S = [histx;histy] + 1e-2*rand(2*n,1);
  		S = [histx;histy] + amplitude*heaviside(t);
	end


set(0,'DefaultAxesFontSize',15)

figure

hAxis(1)=subplot(3,1,1);
pos = get( hAxis(1), 'Position');
pos(1)=.11;
pos(2)=.68;
pos(3)=.87;
pos(4)=.28;
set(hAxis(1), 'Position', pos);
hold on
for ii = 1:n
	plot(sol.x,sol.y(n+ii,:))
end
hold off
box on
ylim([-.4,.4])
set(gca,'YTick',[-.3 0 .3])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gca,'xtick',[]);

hAxis(2)=subplot(3,1,2);
pos = get( hAxis(2), 'Position');
pos(1)=.11;
pos(2)=.37;
pos(3)=.87;
pos(4)=.27;
set(hAxis(2), 'Position', pos);
plot(sol.x,order1,'-b',sol.x,order2,'-r','LineWidth',2);
%ylabel('$R_{1,2}$', txtattrib2{:})
ylim([0,1])
set(gca,'YTick',[0 .5 1])
%ylim([.8,1])
%set(gca,'YTick',[.8 .9 1])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gca,'xtick',[]);

hAxis(3)=subplot(3,1,3);
pos = get( hAxis(3), 'Position');
pos(1)=.11;
pos(2)=.035;
pos(3)=.87;
pos(4)=.28;
set(hAxis(3), 'Position', pos);
hold on
for ii = 1:n
	plot(sol.x(1:end-1),phase_vel(ii,:));
end
hold off
box on
%ylabel('$\hat{\Omega}$', txtattrib2{:})
ylim([-2.5,2.5])
set(gca,'YTick',[-2 0 2])
%ylim([-.12,.12])
%set(gca,'YTick',[-.1 0 .1])
set(gca,'yticklabel',num2str(get(gca,'ytick')','%.1f'))
set(gca,'xtick',[]);

set(gcf, 'PaperPosition', [0 0 8 3]);
set(gcf, 'PaperSize', [8 3]);
saveas(gcf,strcat('trj_',para,'.pdf'));

end
