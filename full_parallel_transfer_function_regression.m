% =========================================================================
% PROFILE FITTING
% Script for fitting real experimental data of microbial fermentation of
% linalool with a parallel effects model. Parallel first order transfer
% functions are employed to approximate the kinetics of terpene
% production upon induction while considering substrate consumption.
% IPTG induction is the system input from 0 to 100% (here equal to 0 to 1). 
% =========================================================================

% =========================================================================
% Initialisation
% =========================================================================
close all
clear all
clc


% =========================================================================
% Data import and preparation
% =========================================================================
% Import data that are saved in a .mat file containing the time values of
% the timepoints ('TimePoints) and the concentration of the target compound
% to be modeled ('Concs'). Both objects are arrays arranged as 1xN, where N
% is the number of the acquired timepoints.
load (['/Users/a30754gl/Desktop/MATLAB works/' ...
    'support scripts/limonene_conc_time_profile.mat'])


% =========================================================================
% Parameters boundaries
% =========================================================================
% Defined over observation of the experimental productivity values
LB = [5 500 500]; % 1st is K, 2nd is Tau1, 3rd is Tau2
UB = [5 4500 4500]; % Same order

iter = 100;
regression_performances = zeros(iter,4);
KoptData = zeros(iter,4);
f = figure('visible','off');
for n = 1:iter
% =========================================================================
% Initialise parameters guesses
% =========================================================================
Kin = LB(1) + (UB(1) - LB(1))/2; 
Tau1in = LB(2) + (UB(2) - LB(2))/2;
Tau2in = LB(3) + (UB(3) - LB(3))/2;


% =========================================================================
% Regress parameters
% =========================================================================
% Minimum search to minimise the objective function, given by Sum of
% Squared Errors (SSE)
X0 = [Kin Tau1in Tau2in];
FOBFun = @(pars)SSECalcFun(pars,Productivity,TimePoints);
options = optimset('MaxIter',100000, ...
    'MaxFunEvals',100000);

[regressed_pars, SSEvalue] = fminsearchbnd(FOBFun,X0, ...
    LB, UB, ...
    options);

% =========================================================================
% Calculate fitted productivity
% =========================================================================
K = regressed_pars(1);
Tau1 = regressed_pars(2);
Tau2 = regressed_pars(3);
DeltaU = 1;
Time = 0:1:4500;

FittedProductivity = ParallelTFStep(K,Tau1,Tau2,DeltaU,Time);


% =========================================================================
% Data plotting
% =========================================================================
% Top graph
if n == 100
    subplot(1,2,1)
    hold on
    plot(Time,FittedProductivity,'b')
    hold on
    scatter(TimePoints,Productivity,'r','x')
    
    legend('Optimal fit','Experimental points')
else
    subplot(1,2,1)
    hold on
    plot(Time,FittedProductivity,'b','HandleVisibility','off')
end

ax = gca;
ax.YLim = [-1 10];


% =========================================================================
% Find K based on current regressed parameters
% =========================================================================
time_index = round(Tau1);
ProdEstimate = FittedProductivity(time_index);
first_effect_contribution = 0.632; % Choosing to pick productivity value at
                                   % 1 time constant (Tau1) means that I am
                                   % at 63.2% of the contribution of the
                                   % first effect by deafult. The
                                   % contribution of the other effect must
                                   % be calculated instead
second_effect_contribution = 1 - exp((-time_index)/(Tau1 + Tau2));

syms X
eqn = ProdEstimate - ...
    (first_effect_contribution*X - second_effect_contribution*X) == 0;
Kestimate = solve(eqn,X);

KoptData(n,:) = ...
    [time_index ProdEstimate second_effect_contribution Kestimate];


% =========================================================================
% Store performance data
% =========================================================================
regression_performances(n,:) = [K Tau1 Tau2 SSEvalue];


% =========================================================================
% Update parameters boundaries
% =========================================================================
UB(1) = UB(1) + 1; % Only K boundaries need to change, as it is the only
                   % parameter which is not dependant on anything else and
                   % can become unphysical
end


% =========================================================================
% Find best fit
% =========================================================================
SSEvalues = regression_performances(:,4);
SSEdiffs = zeros(iter-1);
for i = 2:iter
    SSEdiffs(i-1) = abs(SSEvalues(i-1) - SSEvalues(i));
end

acceptable_indexes = find(SSEdiffs < .01*min(SSEvalues));
best_fit_index = acceptable_indexes(1);


% =========================================================================
% Calculate best profile
% =========================================================================
Kopt = regression_performances(best_fit_index,1);
Tau1opt = regression_performances(best_fit_index,2);
Tau2opt = regression_performances(best_fit_index,3);
FittedProductivity = ParallelTFStep(K,Tau1,Tau2,DeltaU,Time);


% =========================================================================
% Adjust plots
% =========================================================================
% Bottom graph
subplot(1,2,2)
hold on
plot(Time,FittedProductivity,'b')
scatter(TimePoints,Productivity,'r','x')
xlabel('Time [min]')
ylabel('Linalool productivity [mgL^{-1}h^{-1}]')

legend('Optimal fit','Experimental points')

set(f, 'visible', 'on');

disp(regression_performances(best_fit_index,:))


% =========================================================================
% Object function definition (FOBfun)
% ================ =========================================================
% The objective function used in here is a simple Sum of Squared Errors
% (SSE)
function FOBoutput = SSECalcFun(pars,prod,t)
    K = pars(1);
    Tau1 = pars(2);
    Tau2 = pars(3);
    DeltaU = 1;

    y = ParallelTFStep(K,Tau1,Tau2,DeltaU,t);

    FOBoutput = 0;
    for i = 1:length(y)
        FOBoutput = FOBoutput + (prod(i) - y(i))^2;
    end
end


function response_profile = ...
    ParallelTFStep(k,tau1,tau2,delta_u,time)

    response_profile = delta_u*(...
        k*(exp(-time./(tau1 + tau2)) - exp(-time./tau1)));

end









