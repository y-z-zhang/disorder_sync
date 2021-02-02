function stability_code(n,para)
% calculating the maximal transverse Lyapunov exponent (MTLE) of the sync state of delay-coupled Stuart-Landau oscillators by solving characteristic equations
% n -- number of oscillators
% para -- which oscillator to make heterogeneous (choose among 'lambda', 'w', 'gamma')

	% adjacency matrix representing a directed ring network of n nodes
	adj = zeros(1,n);
	adj(1,2) = 1; adj(1,end) = 0;
	ring = adj;
	for ii = 1:n-1
	  adj = [adj;circshift(ring,[0,ii])];
	end
	mu = sum(adj(1,:));
	
	%% parameters for SL oscillators
	beta = 0;
	kmu = .3;
	k = kmu/mu;
	tau = 1.8*pi;
	lambda0 = .1;
	w0 = .12;
	gamma0 = 0;
	h = -.4;
	
	% finding the sync solution
	function y = f_omega(x)
	y = w0 - gamma0*(lambda0 + kmu*(cos(beta-x*tau) - cos(beta)))...
		+ kmu*(sin(beta-x*tau) - sin(beta)) - x;
	end
	
	omega = fzero(@f_omega,w0) % frequency of the sync solution
	phi = beta - omega*tau;
	r2 = lambda0 + kmu*(cos(phi) - cos(beta)) % square amplitude of the sync solution
	

	trial = 10; % number of realizations of random oscillator heterogeneity
	m = 151; sigma = linspace(0,.15,m); % level of heterogeneity measured by standard deviation
	mtle = zeros(trial,m); % initialize maximal transverse Lyapunov exponent
	order = zeros(trial,m); % initialize order parameter

	% main loop
	for jj = 1:trial

		% random seed
		rng('shuffle')
	
		% a random realization of heterogeneity
		sample = normrnd(0,1,n,1);
		sample = (sample - mean(sample))/std(sample);

		% initialize sync state
		r = repmat(sqrt(r2),n,1); Omega = omega; Delta = zeros(n-1,1);
	
		re = .00209;
		% calculate MTLE for each level of heterogeneity
		for ii = 1:m
	
			hetero = sigma(ii); % level of heterogeneity
			dif = sample*hetero; % heterogeneity profile scaled by heterogeneity level

			% adding heterogeneity to the specified parameter(s)
			lambda = lambda0*ones(n,1);
			w = (w0+h)*ones(n,1);
			gamma = (gamma0+h/r2)*ones(n,1);
			switch para
				case 'lambda'
					lambda = lambda + dif;
				case 'w'
					w = w + dif;
				case 'gamma'
					gamma = gamma + dif/r2;
			end
			
			% limit of double-precision numbers ~1e-15
			tol = 1e-15;
			option = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolFun',tol,'TolX',tol);
			options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);
	
			% find the sync state for the given heterogeneous oscillators
			[sol,fval,exitflag,output] = fsolve(@(x)limit_cycle_sol(x,w,gamma,n,k,tau,lambda,kmu,adj),[r;Omega;Delta],options);
			r = sol(1:n); Omega = sol(n+1); delta = [0;sol(n+2:2*n)]; Delta = sol(n+2:2*n);
			if exitflag > 0 % sync state identified successfully
				order(jj,ii) = abs(sum(sign(r).*exp(i*delta)))/n;
			else % no sync state found
				order(jj,ii) = -Inf;
			end
	
			% matrices in the characteristic equation
			J = zeros(2*n);
			for kk = 1:n
				J(2*kk-1:2*kk,2*kk-1:2*kk) = [-2*r(kk)^2 0; -2*gamma(kk)*r(kk)^2 0];
			end
			AR = zeros(2*n);
			P = zeros(2*n);
			for kk = 1:n
				p = zeros(2);
				for kkk = 1:n
					phi = delta(kkk) - delta(kk) - Omega*tau;
					R = r(kkk)/r(kk)*[cos(phi), -sin(phi); sin(phi) cos(phi)];
					AR(2*kk-1:2*kk,2*kkk-1:2*kkk) = adj(kk,kkk)*R;
					p = p + adj(kk,kkk)*R;
				end
				P(2*kk-1:2*kk,2*kk-1:2*kk) = p;
			end
		
			% solving the (transcendental) characteristic equation
			if exitflag > 0 % sync state exist
				floquets = zeros(1,2);
				for kk = 1:1
					[floquets(kk,:),fval,exitflag2,output] = fminsearch(@(x)log_norm_det(x,J,k,P,n,AR,tau),[re,0],option);
				end
				lya_exp = unique(round(floquets(:,1),6))';
				if lya_exp(end) == 0
					floquets = zeros(120,2);
					for kk = 1:120
						[floquets(kk,:),fval,exitflag2,output] = fminsearch(@(x)log_norm_det(x,J,k,P,n,AR,tau),rand(1,2)-.5,option);
					end
					lya_exp = unique(round(floquets(:,1),6))';
					mtle(jj,ii) = lya_exp(end-1);
					re = lya_exp(end-1);
				else
					mtle(jj,ii) = lya_exp(end);
					re = lya_exp(end);
				end 
			else % sync state does not exist
				mtle(jj,ii) = Inf;
			end
		end
	end

	function y = log_norm_det(x,J,k,P,n,AR,tau)
	%% characteristic equation for Floquet exponents
		y = log10(norm(det(J - k*P - (x(1)+x(2)*i)*eye(2*n) + k*AR*exp(-(x(1)+x(2)*i)*tau))));
	end

	function y = limit_cycle_sol(x,w,gamma,n,k,tau,lambda,kmu,adj)
	%% frequency sync solution of nonidentical SL oscillators (works for arbitrary networks with common in-degrees)
		for j = 1:n
			if j == 1
				y(j) = w(j) - gamma(j)*x(j)^2 + k*adj(j,:)*(x(1:n)/x(j).*sin([0;x(n+2:2*n)] - x(n+1)*tau)) - x(n+1);
				y(n+j) = lambda(j) + k*adj(j,:)*(x(1:n)/x(j).*cos([0;x(n+2:2*n)] - x(n+1)*tau)) - kmu - x(j)^2;
			else
				y(j) = w(j) - gamma(j)*x(j)^2 + k*adj(j,:)*(x(1:n)/x(j).*sin([0;x(n+2:2*n)] - x(n+j) - x(n+1)*tau)) - x(n+1);
				y(n+j) = lambda(j) + k*adj(j,:)*(x(1:n)/x(j).*cos([0;x(n+2:2*n)] - x(n+j) - x(n+1)*tau)) - kmu - x(j)^2;
			end
		end
	end

	% save the MTLE and order parameter
	dlmwrite(strcat('mtle_',para,int2str(n),'.dat'),mtle,'-append');
	dlmwrite(strcat('order_',para,int2str(n),'.dat'),order,'-append');

end
