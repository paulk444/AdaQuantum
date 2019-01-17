% Testing script
close all
clear all

%% test making Fock and Coherent states
disp('------------------------------ Basic state making ------------------------------')
%make Fock state |2> truncated at 30
disp('Fock state |2> truncated at 30:')
disp(FockState(2,30))

%make Coherent state alpha = 1 truncated at 3
disp('Coherent state alpha = 1 truncated at 3:')
disp(CoherentState([1,0],3))

%% test photon number probability distribution of coherent state
disp('------------------------------ Coherent states and displacement operator ------------------------------')
alpha = 1 + 2j;
abs_alpha = abs(alpha);
arg_alpha = angle(alpha);
nmax = 25;

probs = zeros(nmax+1,1);
probs2 = zeros(nmax+1,1);
for n = 0:nmax-1
    probs(n+1) = abs(CoherentState([abs_alpha,arg_alpha], nmax)' * FockState(n,nmax))^2;
    probs2(n+1) = abs(alpha)^(2*n)/factorial(n) * exp(- abs(alpha)^2);
end

figure(1)
hold on
scatter(0:nmax, probs)
scatter(0:nmax, probs2,'rx')
%seems to work!
disp('Figure 1 shows photon number probability distribution of coherent states, calculated with |<alpha|i>|^2 and via the formula Barnett + Radmore 3.6.7')
disp(['Checking probability distributions sum to 1. Method 1 gives sum = ' num2str(sum(probs)) '. Method 2 gives sum = ' num2str(sum(probs2))])
fprintf('\n')

%test <alpha|alpha> (will be less than 1 due to truncation error)
disp(['<alpha|alpha> for alpha = ' num2str(alpha) ' and truncation at ' num2str(nmax)])
disp(CoherentState([abs_alpha,arg_alpha], nmax)' * CoherentState([abs_alpha,arg_alpha], nmax))

%% displacement operator

trunc = 15;

%check D(alpha)|0> = |alpha> (will be approx equals due to truncation error)
diff = abs(CoherentState([abs_alpha,arg_alpha], trunc)' * DisplacementOperator([abs(alpha),angle(alpha)], trunc) * FockState(0,trunc));
disp(['Check D(alpha)|0> = |alpha>. For alpha = ' num2str(alpha) ', truncated at ' num2str(trunc) ', <alpha|D(alpha)|0> = ' num2str(diff)])
fprintf('\n')


%% creation and annihilation operators
disp('------------------------------ Creation, annihilation and number operators ------------------------------')

%test a|n> = sqrt(n)|n - 1>. In this case a|3> = sqrt(3)|2> (sqrt(3) =
%1.7321)
disp('a|3> gives:')
disp(AnnihilationOperator(4)*FockState(3,4))

%test a^dag|n> = sqrt(n+1)|n + 1>. In this case a^dag|2> = sqrt(3)|3>
disp('a^dag|2> gives')
disp(CreationOperator(4)*FockState(2,4))

nmax = 20;
%test <alpha|a|alpha> = alpha
disp(['testing <alpha|a|alpha> = alpha. Here alpha = ' num2str(alpha)])
disp(CoherentState([abs_alpha,arg_alpha],nmax)' * AnnihilationOperator(nmax) * CoherentState([abs_alpha,arg_alpha],nmax))

%test <alpha|a^dagger|alpha> = alpha*
disp(['testing <alpha|a^dagger|alpha> = alpha*. Here alpha = ' num2str(alpha)])
disp(CoherentState([abs_alpha,arg_alpha],nmax)' * CreationOperator(nmax) * CoherentState([abs_alpha,arg_alpha],nmax))

%% Number operator

%check <alpha|N|alpha> = |alpha|^2
disp('Testing <alpha|N|alpha> = |alpha|^2.')
disp([' <alpha|N|alpha> = ' num2str(CoherentState([abs_alpha,arg_alpha], trunc)' * NumberOperator(trunc) * CoherentState([abs_alpha,arg_alpha],trunc))]);
disp([' |alpha|^2 = ' num2str(abs(alpha)^2)])
fprintf('\n')

%check <i|N|i> = i
disp('Testing <i|N|i> = i.')
disp(['<5|N|5> = ' num2str(FockState(5, trunc)' * NumberOperator(trunc) * FockState(5, trunc))])
fprintf('\n')

%% Two mode Fock states and number operator

N2 = TMNumberOperator(trunc);
Fock2 = TMFockState(3,4,trunc);

disp(['Testing two mode number operator. <3;4| N | 3;4> = ' num2str(Fock2' * N2 * Fock2)])
fprintf('\n')

%% 1-mode Squeezed states
disp('------------------------------ One mode squeezed states and squeezing operator ------------------------------')
z = 0.5 + 1i;
trunc = 40;
r = abs(z);
phi = angle(z);

%check ~normalised
disp(['Checking normalisation of one mode squeezed state, with z = ' num2str(z) ' truncated at ' num2str(trunc) ':'])
disp(['<z|z> = ' num2str(SingleModeSqueezedState([r,phi], trunc)' * SingleModeSqueezedState([r,phi], trunc))])
fprintf('\n')


probs = zeros(trunc,1);
probs2 = zeros(trunc,1);
for n = 0:trunc-1
    probs(n+1) = abs(SingleModeSqueezedState([r,phi], trunc)' * FockState(n,trunc))^2;
    if mod(n,2) == 0
        m = n/2;
        probs2(n+1) = sech(r) * factorial (2*m)/(factorial(m)^2 * 2^(2*m)) * tanh(r)^(2*m);
    else
        probs2(n+1) = 0;
    end
       
end

figure(2)
hold on
scatter(0:trunc-1, probs,'bo')
scatter(0:trunc-1, probs2,'rx')
disp('Figure 2 shows photon number probability distribution of squeezed state, calculated with |<z|i>|^2 and via the formula Barnett + Radmore 3.7.12-13')
disp(['Checking probability distributions sum to 1. Method 1 gives sum = ' num2str(sum(probs)) '. Method 2 gives sum = ' num2str(sum(probs2))])
fprintf('\n')

disp(['Testing <z|N|z> = sinh(r)^2. For z = ' num2str(z) ', sinh(r)^2 = ' num2str(sinh(r)^2)])
disp(['<z|N|z> = ' num2str(SingleModeSqueezedState([r,phi], trunc)' * NumberOperator(trunc) * SingleModeSqueezedState([r,phi], trunc))])
fprintf('\n')

%% One mode squeezing operator

S = SingleModeSqueezingOperator([r,phi], trunc);
disp(['Testing S(z)|0> = |z>. <z|S(z)|0> = ' num2str(SingleModeSqueezedState([r,phi], trunc)' * S * FockState(0,trunc))])
fprintf('\n')

%% Two mode squeezed state and squeezing operator (p75-76)
disp('------------------------------ Two mode squeezed states and squeezing operator ------------------------------')

disp(['Testing two mode squeezed state and two mode squeezing operator, |z>=S(z)|0>. <z|S(z)|0;0> = '...
num2str(TMSqueezedState([r,phi],trunc)' * TwoModeSqueezingOperator([r,phi],trunc) * TMFockState(0,0,trunc))])
fprintf('\n')

%% Beamsplitter
disp('------------------------------ Beamsplitters ------------------------------')

%action of beamsplitter on Fock states
trunc = 20;

superposition = 1/sqrt(2) * (TMFockState(1,0,trunc) + 1i*TMFockState(0,1,trunc));
superposition2 = 1/sqrt(2) * (TMFockState(1,0,trunc) + TMFockState(0,1,trunc));

disp('Beamsplitter method 1 (book + wiki) has BS|1;0> -> 1/sqrt(2) (|1;0> + i|0;1>). Overlap of these states gives ')
disp(superposition' * BeamSplitter1(0.5, trunc) * TMFockState(1,0,trunc))

disp('Beamsplitter method 2 (Pauls) has BS|1;0> -> 1/sqrt(2) (|1;0> + |0;1>). Overlap of these states gives ')
disp(superposition2' * BeamSplitter(0.5, trunc) * TMFockState(1,0,trunc))

%action of beamsplitter on coherent state
alpha = 0.5 + 0.2i;
abs_alpha = abs(alpha);
arg_alpha = angle(alpha);
beta = 0.1 + 0.3i;
abs_beta = abs(beta);
arg_beta = angle(beta);



T = 0.6;

state1 = BeamSplitter(T, trunc) * TMCoherentState([abs_alpha,arg_alpha], [abs_beta,arg_beta], trunc);
a = sqrt(T)*alpha - sqrt(1 - T)*beta;
b = sqrt(1 - T)*alpha + sqrt(T)*beta;
state2 = TMCoherentState([abs(a),angle(a)], [abs(b),angle(b)] , trunc);
state3 = BeamSplitter1(T, trunc) * TMCoherentState([abs_alpha,arg_alpha], [abs_beta,arg_beta], trunc);
a = sqrt(T)*alpha + 1i*sqrt(1 - T)*beta;
b = 1i*sqrt(1 - T)*alpha + sqrt(T)*beta;
state4 = TMCoherentState([abs(a),angle(a)], [abs(b),angle(b)], trunc);

disp('Beamsplitter method 1 (book + wiki) has BS|alpha; beta> -> |sqrt(T)*alpha + i sqrt(1 - T)*beta ; i sqrt(1 - T)*alpha + sqrt(T)*beta>.')
disp(['For T = ' num2str(T) ', truncated at ' num2str(trunc) '. Overlap of these states gives ' num2str(state4'*state3)])
fprintf('\n')

disp('Beamsplitter method 2 (Pauls) has BS|alpha; beta> -> |sqrt(T)*alpha - sqrt(1 - T)*beta ; sqrt(1 - T)*alpha + sqrt(T)*beta>.')
disp(['For T = ' num2str(T) ', truncated at ' num2str(trunc) '. Overlap of these states gives ' num2str(state2'*state1)])
fprintf('\n')

% simulating annihilation operator with beamsplitter
disp('Simulating annihilation operator with beamsplitter')
epsilon = 0.001;
T = 1-epsilon;
n = 3;

state_produced1 = BeamSplitter1(T,trunc) * TMFockState(n,0,trunc);
state_produced2 = BeamSplitter(T,trunc) * TMFockState(n,0,trunc);
state_predicted1 = TMFockState(n,0,trunc) + 1i*epsilon*sqrt(n)*TMFockState(n-1,1,trunc);
state_predicted2 = TMFockState(n,0,trunc) + epsilon*sqrt(n)*TMFockState(n-1,1,trunc);

disp('Beamsplitter method 1 (book + wiki) has BS|psi; 0> -> |psi;0> + i theta a|psi;1>, for small theta')
disp(['Testing with theta = ' num2str(epsilon) ', truncated at ' num2str(trunc) '. U_BS |' num2str(n) ';0> = |' num2str(n) ';0> + i theta sqrt(' num2str(n) ') |' num2str(n-1) ';1>'])
disp(['< state_produced | state_predicted > = ' num2str(state_produced1' * state_predicted1)])
fprintf('\n')

disp('Beamsplitter method 2 (Pauls) has BS|psi; 0> -> |psi;0> + theta a|psi;1>, for small theta')
disp(['Testing with theta = ' num2str(epsilon) ', truncated at ' num2str(trunc) '. U_BS |' num2str(n) ';0> = |' num2str(n) ';0> + theta sqrt(' num2str(n) ') |' num2str(n-1) ';1>'])
disp(['< state_produced | state_predicted > = ' num2str(state_produced2' * state_predicted2)])
fprintf('\n')

% checking <n> is conserved after beamsplitter for state |alpha,z>
alpha = 0.5 + 0.9i;
z = 0.2;
T = 0.7;

initial_state = kron(CoherentState([abs(alpha),angle(alpha)],trunc),SingleModeSqueezedState([z,0],trunc));
final_state1 = BeamSplitter1(T,trunc) * initial_state;
final_state2 = BeamSplitter(T,trunc) * initial_state;

disp('Checking <n> is conserved before and after beamsplitter. Initial state |alpha,z> should have <n> = |alpha|^2 + sinh^2(r)')
disp(['For alpha = ' num2str(alpha) ', z = ' num2str(z) ', T = ' num2str(T) ', |alpha|^2 + sinh^2(r) = ' num2str(abs(alpha)^2 + sinh(abs(z))^2)])
disp(['<alpha,z|N|alpha,z> = ' num2str(initial_state' * TMNumberOperator(trunc) * initial_state)])
disp(['For method 1, <post-beamsplitter|N|post-beamsplitter> = ' num2str(final_state1' * TMNumberOperator(trunc) * final_state1)])
disp(['For method 2, <post-beamsplitter|N|post-beamsplitter> = ' num2str(final_state2' * TMNumberOperator(trunc) * final_state2)])
fprintf('\n')

%testing TMSS through 50:50 beamsplitter gives SMSS in each mode
trunc = 30;
disp(['Testing TMSS(z) through 50:50 beamsplitter gives |z> in each mode. <z;z| BS | TMSS(z)> = ' num2str(...
kron(SingleModeSqueezedState([z,0],trunc),SingleModeSqueezedState([z,0],trunc))' * BeamSplitter(0.5,trunc) * TMSqueezedState([z,0],trunc)...
)])
fprintf('\n')
%close but not quite - need to find where this is mentioned...

%% Phase shift
disp('------------------------------ Phase Shift ------------------------------')

disp(['Testing phase shift operator. U_phi|1> should give exp(i phi)|1>. <1|U_pi/2|1> should give i. Actually gives: ' num2str(...
FockState(1,trunc)' * PhaseShiftOperator(pi/2,trunc) * FockState(1,trunc)...
)])
fprintf('\n')

%% n^2
disp('------------------------------ n^2 ------------------------------')

alpha = 1 + 0.4i;
trunc = 20;

disp(['Testing <alpha|n^2|alpha> = n(n + 1). <alpha|n^2|alpha> = ' num2str(...
CoherentState([abs(alpha),angle(alpha)],trunc)'*NSquared(trunc)*CoherentState([abs(alpha),angle(alpha)],trunc)...
) '. n(n + 1) = ' num2str(...
abs(alpha)^2*(1 + abs(alpha)^2)...
)])
fprintf('\n')

disp(['Testing <4|n^2|4> = 4^2. <4|n^2|4> = ' num2str(...
FockState(4,trunc)'*NSquared(trunc)*FockState(4,trunc)...
)])
fprintf('\n')

%% QFI
disp('------------------------------ QFI ------------------------------')
trunc = 100;

state = CoherentState([abs(alpha),angle(alpha)],trunc);
disp(['Testing QFI of Coherent state using variance of n method. QFI(|alpha>) = '...
num2str(PureStateQFI(state * state'))...
'. Should = 4n = '...
num2str(4 * abs(alpha)^2)])

state = FockState(3,trunc);
disp(['Testing QFI of Fock state using variance of n method. QFI(|3>) = '...
num2str(PureStateQFI(state * state'))...
'. Should = 0'])
fprintf('\n')

state = CoherentState([abs(alpha),angle(alpha)],trunc);
disp(['Testing QFI of Coherent state using eigenvalues + eigenvectors method. Truncated at '...
num2str(trunc)...    
'. QFI(|alpha>) = '...
num2str(MixedStateQFI(state * state'))...
'. Should = 4n = '...
num2str(4 * abs(alpha)^2)])

state = CoherentState([abs(alpha),angle(alpha)],trunc);
state = 98/100*(state * state') + 2/100*eye(trunc);
disp(['Testing QFI of slightly mixed Coherent state using eigenvalues + eigenvectors method. Truncated at '...
num2str(trunc)...    
'. QFI(rho) = '...
num2str(MixedStateQFI(state))...
'. Should be slightly less than 4n = '...
num2str(4 * abs(alpha)^2)])

state = FockState(3,trunc);
disp(['Testing QFI of Fock state using eigenvalues + eigenvectors method. QFI(|3>) = '...
num2str(MixedStateQFI(state * state'))...
'. Should = 0'])
fprintf('\n')

state = CoherentState([abs(alpha),angle(alpha)],trunc);
rho = state * state';
disp('Testing QFI of Coherent state using general QFI function.')
TheQFI = QFI(rho);
disp(['QFI(|alpha>) = ' num2str(TheQFI) '. Should = 4n = '...
num2str(4 * abs(alpha)^2)])
fprintf('\n')

state = CoherentState([abs(alpha),angle(alpha)],trunc);
rho = 98/100*(state * state') + 2/100*eye(trunc);
disp('Testing QFI of slightly mixed Coherent state using general QFI function.')
TheQFI = QFI(rho);
disp(['QFI(rho) = ' num2str(TheQFI) '. Should be slightly less than 4n = '...
num2str(4 * abs(alpha)^2)])
fprintf('\n')

%% Quadrature operator and eigenstate
disp('------------------------------ Quadrature operator and eigenstate ------------------------------')
trunc = 30;

alpha = 1+2i;
x_lambda = 1 + 1i*sqrt(3);
r = abs(x_lambda);
lambda = angle(x_lambda);

comm = Commutator(QuaderatureOperator(lambda, trunc), QuaderatureOperator(lambda+pi/2,trunc));
disp('Testing commutator of x_lambda and x_lambda+pi/2 = i*Identity. Have to discount last element due to truncation error. Average diagonal element = ')
disp(trace(comm(1:trunc-1,1:trunc-1))/(trunc-1))

%disp(['Testing both methods of creating quadrature eigenstates. Sum of difference between the two = '...
%num2str(sum(abs(QuaderatureEigenstate([r,lambda], trunc) - QuaderatureEigenstate2(x_lambda, trunc))))])
%fprintf('\n')

disp('Testing eigenvalue equation. Have to discount last element due to truncation error')
should_be = QuaderatureOperator(lambda, trunc+1)*QuaderatureEigenstate([r,lambda], trunc+1);
actually = r * QuaderatureEigenstate([r,lambda], trunc+1);
total_error = sum(abs(should_be(1:trunc) - actually(1:trunc)));
disp(['Total of errrors x^hat_lambda |x_lambda> - x_lambda |x_lambda> = ' num2str(total_error)])
fprintf('\n')

disp(['Testing <x_lambda|alpha> against A.4.12 in the book. <x_lambda|alpha> gives '...
num2str(QuaderatureEigenstate([r,lambda], trunc)' * CoherentState([abs(alpha),angle(alpha)],trunc))...
'. Should be '...
num2str(QuadratureCoherentOverlap(x_lambda, alpha))])
fprintf('\n')

trunc = 10;
dx = 0.1;
lambda = pi/4;
sum = zeros(trunc);
for x_lambda = -100:dx:100
    quadrature_eig = QuaderatureEigenstate([x_lambda,lambda], trunc);
    projection_operator = quadrature_eig * quadrature_eig';
    sum = sum + projection_operator * dx;
end
disp('Checking Sum( |x_lambda><x_lambda| dx) = Identity. Diagonal of this sum = ')
disp(diag(sum))

% trunc_max = 150;
% Method1times = zeros(trunc_max,1);
% Method2times = zeros(trunc_max,1);
% 
% for trunc = 1:trunc_max
%     f1 = @() QuaderatureEigenstate(x_lambda, trunc);
%     f2 = @() QuaderatureEigenstate2(x_lambda, trunc);
%     Method1times(trunc) = timeit(f1);
%     Method2times(trunc) = timeit(f2);
% end
% 
% figure
% hold on
% scatter(1:trunc, Method1times,'bo')
% scatter(1:trunc, Method2times,'rx')
disp('Tested timing Method 1 vs. Method 2. Method 1 is MUCH quicker')
fprintf('\n')
%% -------------------- MEASUREMENTS --------------------------

%% Number measurement
disp('------------------------------ Number and Bucket Measurements ------------------------------')
trunc = 5;

state_in = (3*TMFockState(0,2,trunc) + TMFockState(0,1,trunc) + 2*TMFockState(1,2,trunc) + 4*TMFockState(2,1,trunc))/sqrt(30);

for i = 0:2
[output_state, output_probability] = NumberMeasurement(state_in, i);
alternate_output_state = NumberMeasurementToCompare(state_in, i);

disp(['When doing Fock measurement |' num2str(i) '><' num2str(i) '| Output state = '])
disp(output_state)
disp(['with probability, p = ' num2str(output_probability)])
disp('Using alternate method, outputstate = ')
%disp(output_state_test);
%disp(['with probability, p = ' num2str(output_probability_test)])

end
[output_state, output_probability] = BucketMeasurement(state_in, true);
disp('When doing Bucket measurement. Output state = ')
disp(output_state)
disp(['with probability, p = ' num2str(output_probability)])

disp('All results as they should be')
fprintf('\n')

%% Homodyne measurement
%testing in ReproduceMolnar2017.m

%% Density matrix to state vector
alpha = 0.6 + 0.8i;
trunc = 30;

state = CoherentState([abs(alpha),angle(alpha)], trunc);
density_matrix = state * state';
state2 = DensityOperatorToStateVector(density_matrix);
overlap = state2' * state;

disp(['Testing DensityOperatorToStateVector reproduces input coherent state. Overlap between original and output = ' num2str(overlap)])
fprintf('\n')

%% Measurements on arbitrary numbers of modes



