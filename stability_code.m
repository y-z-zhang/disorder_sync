function mtle(para,n)
%% calculate the MTLE of a directed ring of n delay-coupled
%% Stuart-Landau oscillators analytically

	% adjaency matrix representing ring netwrok
	adj = zeros(1,n);
	adj(1,2) = 1; adj(1,end) = 0;
	ring = adj;
	for ii = 1:n-1
	  adj = [adj;circshift(ring,[0,ii])];
	end
	mu = sum(adj(1,:));

	% adjaency matrix representing regular random netwrok
	%adj = dlmread('reg_ran_2.txt');
	%mu = sum(adj(1,:));

	% adjaency matrix representing an open 2D lattice
	%mu = 4;
	%ring = zeros(1,n);
	%ring(1,2) = 1; ring(1,end) = 0;
	%adj = toeplitz(ring);
	%adj = kron(eye(n),adj) + kron(adj,eye(n));
	%for ii = 1:n^2
	%  adj(ii,:) = mu/sum(adj(ii,:))*adj(ii,:);
	%end
	%n = n^2;
	
	%% parameters for SL oscillators
	beta = 0;
	kmu = .3; 
	%kmu = .1;
	k = kmu/mu;
	%tau = .8*pi;
	tau = 1.8*pi;
	lambda0 = .1;
	%w0 = .1;
	w0 = .12;
	%w0 = 1;
	gamma0 = 0;
	
	function y = f_omega(x)
	%% solve omega as x
	y = w0 - gamma0*(lambda0 + kmu*(cos(beta-x*tau) - cos(beta)))...
		+ kmu*(sin(beta-x*tau) - sin(beta)) - x;
	end
	
	omega = fzero(@f_omega,w0)
	
	phi = beta - omega*tau;
	
	r2 = lambda0 + kmu*(cos(phi) - cos(beta))
	
	h = -.4;
	trial = 1000;
	m = 301;
	heterogeneity = linspace(0,.3,m);
	mtle = zeros(trial,m);
	order = zeros(trial,m);

	parfor jj = 1:trial

		rng('shuffle')
		%rng(jj)
	
		sample = normrnd(0,1,n,1);
		sample = (sample - mean(sample))/std(sample);
		sample1 = normrnd(0,1,n,1);
		sample1 = (sample1 - mean(sample1))/std(sample1);
		sample2 = normrnd(0,1,n,1);
		sample2 = (sample2 - mean(sample2))/std(sample2);

		r = repmat(sqrt(r2),n,1); Omega = omega; Delta = zeros(n-1,1);
	
		re = .00209;
		for ii = 1:m
	
			hetero = heterogeneity(ii);	
			dif = sample*hetero;
			dif1 = sample1*hetero;
			dif2 = sample2*hetero;
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
				case 'special'
					w = w + dif;
					gamma = gamma + dif/r2;
				case 'w+gamma'
					w = w + dif;
					gamma = gamma + dif1/r2;
				case 'all'
					lambda = lambda + dif;
					w = w + dif1;
					gamma = gamma + dif2/r2;
			end
			
			%% limit of double-precision numbers ~1e-15
			tol = 1e-15;
			option = optimset('Display','off','MaxIter',1000,'MaxFunEvals',1000,'TolFun',tol,'TolX',tol);
			options = optimoptions('fsolve','Display','off','MaxIterations',1000,'MaxFunctionEvaluations',100000);%,...
				%'FunctionTolerance',tol,'StepTolerance',tol,'OptimalityTolerance',tol);
	
			[sol,fval,exitflag,output] = fsolve(@(x)limit_cycle_sol(x,w,gamma,n,k,tau,lambda,kmu,adj),[r;Omega;Delta],options);
			r = sol(1:n); Omega = sol(n+1); delta = [0;sol(n+2:2*n)]; Delta = sol(n+2:2*n);
			if exitflag > 0
				order(jj,ii) = abs(sum(sign(r).*exp(i*delta)))/n;
			else
				order(jj,ii) = -Inf;
			end
	
		
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
		
			if exitflag > 0
				floquets = zeros(1,2);
				for kk = 1:1
					%[floquets(kk,:),fval,exitflag2,output] = fsolve(@(x)characteristic_eq(x,J,k,P,n,AR,tau),[re,0],options);
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
			else
				mtle(jj,ii) = Inf;
			end
		end
	end

	%% bootstrapping
	%success = (mtle(ii,:)<0).*(order(ii,:)>.7);
	%successes = bootstrp(mm,@mean,success);
	%ave = mean(success)
	%err = std(successes)
	%improve = (mtle(ii,:)<mtle(1,1)).*(order(ii,:)>.7);
	%improves = bootstrp(mm,@mean,improve);
	%ave1 = mean(improve)
	%err1 = std(improves)

	dlmwrite(strcat('mtle_',para,int2str(n),'.dat'),mtle,'-append');
	dlmwrite(strcat('order_',para,int2str(n),'.dat'),order,'-append');

end
