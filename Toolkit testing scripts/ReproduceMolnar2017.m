function ReproduceMolnar2017(line)
% true for line, false for lattice
% Testing toolbox against arXiv:1708.02896v1 2017, Molnar et al.
% "Quantum state engineering via coherent-state superpositions in traveling optical fields"
startup

%% Choose truncation and tolerance
trunc = 40;
epsilon_normalisation = 0.01;

%% Input states
phi = rand() * pi/2;

abs_alpha = rand() * 2;
arg_alpha = zeros(2,1);
arg_alpha(1) = pi/2 + phi/2;
if line
    arg_alpha(2) = pi/2 + phi/2;
else
    arg_alpha(2) = phi/2;
end
alpha = abs_alpha * exp(1i*arg_alpha);

normalisation_const = zeros(2,1);
psi_in = zeros(trunc,2);
input = zeros(trunc^2, 2);
for i = 1:2
    normalisation_const(i) = 1/sqrt(2 + exp(-abs(alpha(i))^2) * (exp(abs(alpha(i))^2 * exp(-1i*phi)) + exp(abs(alpha(i))^2 * exp(-1i*phi))));
    psi_in(:,i) = normalisation_const(i) * (CoherentState(alpha(i),trunc) + CoherentState(alpha(i)*exp(-1i*phi),trunc));
    overlap = psi_in(:,i)' * psi_in(:,i);
    disp(['Overlap for psi_in(' num2str(i) ') = ' num2str(overlap)])
    if abs(1 - overlap) > epsilon_normalisation
        disp(['WARNING: State not normalised. Overlap for psi_in(' num2str(i) ') = ' num2str(overlap)])
    end
    input(:,i) = kron(psi_in(:,i),psi_in(:,i));
end


%% Define gates

BS = BeamSplitter(0.5, trunc);

homodyne_mag = zeros(3,1);
for i = 1:3
    homodyne_mag(i) = rand() * 1;
end

if line
    homodyne_phases = [0;0;0];
else
    homodyne_phases = [0;pi/2;0];
end
homodyne_settings = homodyne_mag .* exp(1i * homodyne_phases);


%% Act gates

FirstHomodyneOut = zeros(trunc,trunc,2);
for i = 1:2
    FirstHomodyneOut(:,i) = DensityOperatorToStateVector(HomodyneMeasurementOnMode2( BS * input(:,i), homodyne_settings(i)));
end


CombinedFirstHomodyneOut = kron(FirstHomodyneOut(:,1),FirstHomodyneOut(:,2));

OutputState = DensityOperatorToStateVector( HomodyneMeasurementOnMode2( BS * CombinedFirstHomodyneOut, homodyne_settings(3) ) );

%% Mid-point theoretical output

rt2alpha = sqrt(2)*abs_alpha;

% Eq. 8
a_coeffs = zeros(2,1);
a_coeffs(1) = QuadratureCoherentOverlap(homodyne_settings(1), rt2alpha*1i*exp(1i*phi/2)) +  QuadratureCoherentOverlap(homodyne_settings(1), rt2alpha*1i*exp(-1i*phi/2)); %a_0
a_coeffs(2) = QuadratureCoherentOverlap(homodyne_settings(1), rt2alpha*1i*cos(phi/2)); %a_1
b_coeffs = zeros(2,1);
b_coeffs(1) = QuadratureCoherentOverlap(homodyne_settings(2), rt2alpha*1i*exp(1i*phi/2)) + QuadratureCoherentOverlap(homodyne_settings(2), rt2alpha*1i*exp(-1i*phi/2)) ; %b_0
b_coeffs(2) = QuadratureCoherentOverlap(homodyne_settings(2), rt2alpha*1i*cos(phi/2)); %b_1

% Eq. 7
cat_state = CoherentState(rt2alpha*sin(phi/2),trunc) + CoherentState( - rt2alpha*sin(phi/2),trunc);
cat_prime_state = CoherentState(rt2alpha*1i*sin(phi/2),trunc) + CoherentState( - rt2alpha*1i*sin(phi/2),trunc);

psi_1_mid = a_coeffs(1) * FockState(0, trunc) + a_coeffs(2) * cat_state;
if line
    psi_2_mid = b_coeffs(1) * FockState(0,trunc) + b_coeffs(2) * cat_state;
else
    psi_2_mid = b_coeffs(1) * FockState(0,trunc) + b_coeffs(2) * cat_prime_state;
end

% normalise
overlap1 = psi_1_mid' * psi_1_mid;
overlap2 = psi_2_mid' * psi_2_mid;
normalisation_const_1 = 1/sqrt(overlap1);
normalisation_const_2 = 1/sqrt(overlap2);
theoretical_mid_state_1 = normalisation_const_1 * psi_1_mid;
theoretical_mid_state_2 = normalisation_const_2 * psi_2_mid;

%% Compare mid-point output to theory

disp(['Overlap between theoretical mid state (1) and mid state (1) through gates = ' num2str(abs(theoretical_mid_state_1' * FirstHomodyneOut(:,1)))])

figure
plot(0:trunc-1,abs(FirstHomodyneOut(:,1)), 'rx')
hold on
plot(0:trunc-1,abs(theoretical_mid_state_1), 'bo')
hold off

disp(['Overlap between theoretical mid state (2) and mid state (2) through gates = ' num2str(abs(theoretical_mid_state_2' * FirstHomodyneOut(:,2)))])

figure
plot(0:trunc-1,abs(FirstHomodyneOut(:,2)), 'rx')
hold on
plot(0:trunc-1,abs(theoretical_mid_state_2), 'bo')
hold off

%% Theoretical output

beta = abs_alpha * sin(phi/2);

theoretical_output_state = zeros(trunc,1);
if line
   
    % Eq. 9
    coeffs = zeros(5,1);
    coeffs(1) = a_coeffs(2) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), 0); %c_-2
    coeffs(2) = a_coeffs(1) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta) + a_coeffs(2) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), -beta); %c_-1
    coeffs(3) = a_coeffs(1) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), 0) + a_coeffs(2) * b_coeffs(2) * (QuadratureCoherentOverlap(homodyne_settings(3), 2*beta) + QuadratureCoherentOverlap(homodyne_settings(3), -2*beta)); %c_0
    coeffs(4) = a_coeffs(1) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), -beta) + a_coeffs(2) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), beta); %c_1
    coeffs(5) = coeffs(1); %c_2
    
    for n = -2:2
        theoretical_output_state = theoretical_output_state + coeffs(n+3) * CoherentState( n * beta, trunc);
    end
    
    overlap = theoretical_output_state' * theoretical_output_state;
    normalisation_const = 1/sqrt(overlap);
    theoretical_output_state = normalisation_const * theoretical_output_state;
    
else
    
    % Eq. 10
    coeffs = zeros(3);
    coeffs(1,3) = a_coeffs(2) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta*(-1-i)); %c_-1,1
    coeffs(2,3) = a_coeffs(1) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), -beta); %c_0,1
    coeffs(3,3) = a_coeffs(2) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta*(1-i)); %c_1,1
    coeffs(1,2) = a_coeffs(2) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), -beta); %c_-1,0
    coeffs(2,2) = a_coeffs(1) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), beta); %c_0,0
    coeffs(3,2) = a_coeffs(2) * b_coeffs(1) * QuadratureCoherentOverlap(homodyne_settings(3), beta); %c_1,0
    coeffs(1,1) = a_coeffs(2) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta*(-1+i)); %c_-1,-1
    coeffs(2,1) = a_coeffs(1) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta*1i); %c_0,-1
    coeffs(3,1) = a_coeffs(2) * b_coeffs(2) * QuadratureCoherentOverlap(homodyne_settings(3), beta*(1+i)); %c_1,-1
    
    
    for k = -1:1
        for l = -1:1
            theoretical_output_state = theoretical_output_state + coeffs(k+2,l+2) * CoherentState( (k + 1i*l) * beta, trunc);
        end
    end
    
    overlap = theoretical_output_state' * theoretical_output_state;
    normalisation_const = 1/sqrt(overlap);
    theoretical_output_state = normalisation_const * theoretical_output_state;
    
end

%% Compare result to theory

disp(['Overlap between theoretical output state and output state through gates = ' num2str(abs(theoretical_output_state' * OutputState))])

figure
plot(0:trunc-1,abs(OutputState), 'rx')
hold on
plot(0:trunc-1,abs(theoretical_output_state), 'bo')




