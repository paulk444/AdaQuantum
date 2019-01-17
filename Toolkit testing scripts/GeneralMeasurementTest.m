% Script to test the measurement functions (after introduction of general
% measurement function)
close all
clear all

%% Two-mode input, one-mode output - comparison with old measurement functions

%% Number measurement
disp('------------------------------ Number Measurements ------------------------------')
trunc = 5;

state_in = (3*TMFockState(0,2,trunc) + TMFockState(0,1,trunc) + 2*TMFockState(1,2,trunc) + 4*TMFockState(2,1,trunc))/sqrt(30);
rho_in = state_in * state_in';

for i = 0:2
[output_state, output_probability] = NumberMeasurement(state_in, i);
[output_state_new, output_probability_new] = GeneralMeasurement(rho_in, {NumberPOVM(i, 0, trunc)});

disp(['When doing Fock measurement |' num2str(i) '><' num2str(i) '| Output state = '])
disp(output_state)
disp(['with probability, p = ' num2str(output_probability)])
disp('Using new method, outputstate = ')
disp(output_state_new);
disp(['with probability, p = ' num2str(output_probability_new)])

end

%% Bucket measurements

[output_state, output_probability] = BucketMeasurement(state_in, true);
[output_state_new, output_probability_new] = GeneralMeasurement(rho_in, {BucketPOVM(true, 0, trunc)});

disp('When doing Bucket measurement. Output state = ')
disp(output_state)
disp(['with probability, p = ' num2str(output_probability)])

disp('Using new method, outputstate = ')
disp(output_state_new);
disp(['with probability, p = ' num2str(output_probability_new)])

%% Homodyne?

%% Multiplex detector

trunc = 10;

%Input = kron(FockState(2,trunc), CoherentState([1,0],trunc)); 
Input = NormaliseState(TMFockState(1,1,trunc) + TMFockState(2,2,trunc));
Rho = Input * Input';

[Output, Prob] = MultiplexDetector(Input, 4, 1);

[OutputG, ProbG] = GeneralMeasurement(Rho, {MultiplexDetectorPOVM(1,[4,0],trunc)});
[OutputP, ProbP] = POVMMeasurement(Input, MultiplexDetectorPOVM(1,[4,0],trunc));

{'Prob', Prob; 'ProbG', ProbG; 'ProbP', ProbP}

%% Compare timings projectors vs. lossy POVMs

trunc = 35;

state_in = (3*TMSqueezedState([1.4,1.11],trunc) + TMFockState(0,1,trunc) + 2*TMCoherentState([1,0],[2,1.32],trunc) + 4*TMFockState(2,1,trunc))/sqrt(30);
rho_in = state_in * state_in';

projTimes = zeros(3,1);
povmTimes = zeros(3,1);
for i = 0:2
    tic
    [output_state_proj, output_probability_proj] = GeneralMeasurement(rho_in, {NumberProjector(i, trunc)});
    projTimes(i+1) = toc;
    tic
    [output_state_new, output_probability_new] = GeneralMeasurement(rho_in, {NumberPOVM(i, 0, trunc)});
    povmTimes(i+1) = toc;
end

tic
[output_state_new, output_probability_new] = GeneralMeasurement(rho_in, {BucketProjector(true, trunc)});
bucketproj = toc;
tic
[output_state_new, output_probability_new] = GeneralMeasurement(rho_in, {BucketPOVM(true, 0, trunc)});
bucketpovm = toc;

% They seem to be roughly comparable

%% compare paul's vs mine speed

%% Compare timings projectors vs. lossy POVMs

trunc = 80;

state_in = (3*TMSqueezedState([1.4,1.11],trunc) + TMFockState(0,1,trunc) + 2*TMCoherentState([1,0],[2,1.32],trunc) + 4*TMFockState(2,1,trunc))/sqrt(30);
rho_in = state_in * state_in';

paulTimes = zeros(3,1);
rosannaTimes = zeros(3,1);
for i = 0:2
    tic
    [output_state_proj, output_probability_proj] = GeneralMeasurement(rho_in, {NumberPOVM(i, 0, trunc)});
    paulTimes(i+1) = toc;
    tic
    [output_state_new, output_probability_new] = GeneralMeasurement_PKV(rho_in, {NumberPOVM(i, 0, trunc)});
    rosannaTimes(i+1) = toc;
end
