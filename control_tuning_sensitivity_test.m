% =========================================================================
% PID TUNING SENSITIVITY TEST
% PID controller is tuned in order to execute feedback control on linalool
% productivity. Biological system transfer function is derived from
% experimental data and used here to build a close loop control system.
% Control efficacy sensitivity is tested by modifying PID parameters up to
% 20% from their optimal value in a random way. 1000 scenarios are randomly
% generated.
% =========================================================================

% =========================================================================
% Initialisation
% =========================================================================
close all
clear all
clc

format long
rng('shuffle');


% =========================================================================
% Closed loop building
% =========================================================================
% Definition of bioprocess transfer function from previously obtained
% parameters
numerator = 23;
denominator = [1219.6 1];
sys = tf(numerator,denominator);
sys.TimeUnit = 'minutes';

% PID is tuned on the basis of the bioprocess transfer function
[C_pid,info] = pidtune(sys,'PIDF',0.5);
C_pid.TimeUnit = 'minutes';

% Outer cycle defining 2 cases: step change & disturbance rejection
for c = 1:2
    if c == 1
        CL_pid = feedback(C_pid*sys,1);
        [y1,t1] = step(CL_pid);
        
        figure(c)
        hold on
    else
        CL_pid_dist = feedback(sys,C_pid);
        [y1,t1] = step(CL_pid_dist);
        
        figure(c)
        hold on
    end

    % Inner cycle to generate 1000 random combinations of the 4 parameters,
    % varied up to a maximum of 20% from thei optimal value
    pid_tunings = zeros(1001,4);
    pid_tunings(1,:) = [C_pid.Kp C_pid.Ki C_pid.Kd C_pid.Tf];
    ranges = abs(.20*pid_tunings(1,:));
    for i = 1:1000
        for n = 1:4
            eval(['random' num2str(n) ' = rand;']);
            eval([ ...
            'if random' num2str(n) ' >=.5;' ...
                 'sign' num2str(n) ' = 1;' ...
            'else;' ...
                 'sign' num2str(n) ' = -1;' ...
            'end']);
        end

        C_pid.Kp = pid_tunings(1,1) + sign1*random1*ranges(1);
        C_pid.Ki = pid_tunings(1,2) + sign2*random2*ranges(2);
        C_pid.Kd = pid_tunings(1,3) + sign3*random3*ranges(3);
        C_pid.Tf = pid_tunings(1,4) + sign4*random4*ranges(4);

        pid_tunings(n+1,:) = [C_pid.Kp C_pid.Ki C_pid.Kd C_pid.Tf];
        
        if c == 1
            CL_pid = feedback(C_pid*sys,1);
        else
            CL_pid = feedback(sys,C_pid);
        end

        [y,t] = step(CL_pid);
        if i ~= 1000
            plot(t,y,'r','HandleVisibility','off')
        else
            plot(t,y,'r')
        end

    end

    plot(t1,y1,'b')
    xlabel('Minutes')
    ylabel('\DeltaProductivity/\DeltaProductivity_{max}')
    legend('Random tuning parameters deviation (up to 20% from optimal)', ...
        'Optimal tuning parameters')

end
