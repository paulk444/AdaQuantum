function [bayes_bound] = mz_1trial_bound_evo(initial_state,phase_width,phase_mean)
% Jesús Rubio Jiménez, PhD student
% University of Sussex
% 27th August 2018
% J.Rubio-Jimenez@sussex.ac.uk
%
% mz_1trial_bound_evo(initial_state,phase_width,phase_mean),
%
% where 'initial_state' is the initial state of a Mach-Zehnder interferometer,
% 'phase_width' is the width of a flat prior density function and 'phase_mean'
% is the mean of the same flat prior.
%
% Following the quantum Bayesian techniques developed by Helstrom and 
% Personick, this programme calculates the value of the optimal Bayesian
% mean square error 'bayes_bound' for 1 observation/trial/repetition.
%
% Note that mz_optimal_1trial_evo also calculated 
%
%   - the optimal Hermitian observable 'Lopt',
%   - the possible estimates 'lambda' for the unknown parameter given by
%     the spectrum of Lopt,
%   - the optimal projective measurement for 1 observation/trial/repetition 
%     given by the eigenvectors 'lambdavec_columns' of Lopt,
%   - the zero-th moment of the transformed density matrix multiplied by
%     the prior probability 'rho',
%   - the diagonal zero-th moment of the transformed density matrix multiplied
%     by the prior probability 'pk', 
%   - the matrix 'psik' whose columns are the eigenvectors of rho,
%     the first moment of the transformed density matrix multiplied by
%     the prior probability 'rhobar',
%   - the first moment of the transformed density matrix multiplied by
%     the prior probability 'rhobarnew' expressed in the eigenbasis of rho.
%
% These outputs are removed in order to increase the speed for the calculation 
% of the optimal error itself, which is the only output that we need. 

% ADDED by Jesús Rubio on Jan 2021: for the theory behind the optimal quantum  
% mean square error as calculated here, see
%   * Rubio, J., and Dunningham, J. (2019). Quantum metrology in the presence
%     of limited data. New Journal of Physics, 21 043037
%   * Rubio Jiménez, J. (2020). Non-asymptotic quantum metrology: extracting
%     maximum information from limited data. PhD thesis, University of Sussex,
%     ISNI: 0000 0004 8504 6357 (arXiv:1912.02324).

% Bayes statistical and mean operators
index=1;
kvec=zeros(1,length(initial_state)^2);
lvec=zeros(1,length(initial_state)^2);
for x1=1:sqrt(length(initial_state))
    for y1=1:sqrt(length(initial_state))
        for z1=1:sqrt(length(initial_state))
            for t1=1:sqrt(length(initial_state))
                if (x1-1)-(y1-1)+(t1-1)-(z1-1)==0
                    K=phase_width;
                    L=phase_mean*phase_width;
                else
                    comp_temp=(x1-1)-(y1-1)+(t1-1)-(z1-1);
                    exp_temp=exp(-1i*comp_temp*phase_mean/2);
                    sin_temp=sin(comp_temp*phase_width/4);
                    cos_temp=cos(comp_temp*phase_width/4);
                    K=4*exp_temp*sin_temp/comp_temp;
                    L=exp_temp*(4*phase_mean*sin_temp/comp_temp + 1i*2*phase_width*cos_temp/comp_temp - 1i*8*sin_temp/comp_temp^2);
                end
                kvec(index)=K/phase_width;
                lvec(index)=L/phase_width;
                index=index+1;
            end
        end
    end
end
kmat=sparse(vec2mat(kvec,sqrt(length(kvec))));
lmat=sparse(vec2mat(lvec,sqrt(length(lvec))));
initial_rho=kron(initial_state,initial_state');
rho=initial_rho.*kmat;
rhobar=initial_rho.*lmat;
rho=full(rho);
rhobar=full(rhobar);

% Eigenvalues and eigenvectors for the Bayes statistical operator
[psik, pk] = eigs(rho,rank(rho));
psik=sparse(psik);
pk=sparse(pk);

pkvec=zeros(1,length(pk));
for x=1:length(pk)
    pkvec(x)=pk(x,x);
end

% Bayes mean operator in the basis of eigenvectors of the Bayes statistical operator
rhobarnew=psik'*rhobar*psik;

% Optimal projective measurements and its experimental outcomes
Lopt_temp=zeros(length(pkvec),length(pkvec));
for a=1:length(pkvec)
    for b=1:length(pkvec)
        if pkvec(a)+ pkvec(b)>0
            Lopt_temp(a,b)=2*rhobarnew(a,b)/(pkvec(a)+pkvec(b));
        end
    end
end
Lopt_temp=sparse(Lopt_temp);
Lopt=psik*Lopt_temp*psik';
Lopt=full(Lopt);

% Bayes optimal bound for 1 observation/trial/repetition
bayes_bound=phase_width^2/12+phase_mean^2-trace(Lopt*Lopt*rho);
if imag(bayes_bound)<1e-10
    bayes_bound=real(bayes_bound);
else
    error('The mean square error must be a real quantity.')
    return
end
%disp('The optimal POVM for 1 observation has been created.')
end