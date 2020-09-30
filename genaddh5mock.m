function name_curr_mock = genaddh5mock(user_mass,timelength_mock,varargin)
% =========================================================================
% INPUTs
% 'user_mass' = Mass value at which simulate a random intensity profile
% 'timelength_mock' = Timespan covered by the mockfile to be generated
% 
% OPTIONAL INPUTs
% (1st ARG)'random distribution function' = The function used to generate
%                                           the peak intensity profile over
%                                           time.
%                                           DEFAULT = 'stable'
%                                           OPT1 = 'sinusoindal'
%                                           OPT2 = 'constant'
% (2nd ARG)'measured masses starting point' = The lowest mass value
%                                             detected by the PTRMS during
%                                             the current mock.
%                                             DEFAULT = 0
%                                             OPT = N > 0
% 
% OUTPUTs
% 'name_curr_mock' = Name of the mock .h5 file in output
% 
%
% Function that generates a mock .h5 file to emulate the output of the
% PTR-MS. The output file is generated exactly as the original one, except
% for the data logging regarding the number of completed acquisition cycles
% and the number of zeroes contained in the last time cycle. In addition to
% file generation, the function adds a random output profile at the mass
% value chosen by the user. This random profile is given by the combination
% of two distributions: (i) a distribution over time based on a userdefined
% distribution' which determines the intensity profile over time for the
% exact mass inputted by the user; (ii) a normal distribution across the
% masses to generate fictional peaks around the chosen mass (usually there
% is a neighbourhood of masses to be considered when detecting a specific
% one)
% =========================================================================


% =========================================================================
% Initialisation
% =========================================================================
% delete 'PTRMSmocksequence*.h5'
% close all
% clear all
format long
% clc


% =========================================================================
% Error handling & optional arguments evaluation
% =========================================================================
switch nargin
    case 2
        assert(isnumeric(user_mass) & ...
            2 == ndims(user_mass) & ...
            216 > user_mass & ...
            0 < user_mass & ...
            1 == size(user_mass,1) & ...
            1 == size(user_mass,2), ...
            ['First input <user_mass> must be a positive scalar ' ...
            'lower than 216m/q.'])
        assert(isnumeric(timelength_mock) & ...
            2 == ndims(timelength_mock) & ...
            5 < timelength_mock & ...
            1 == size(timelength_mock,1) & ...
            1 == size(timelength_mock,2), ...
            ['Second input <timelength_mock> must be a positive ' ...
            'scalar higher than 5s.'])
        
        start_mass = 0;
        distr_fun = 'stable';
    case 3
        assert(isnumeric(user_mass) & ...
            2 == ndims(user_mass) & ...
            216 > user_mass & ...
            0 < user_mass & ...
            1 == size(user_mass,1) & ...
            1 == size(user_mass,2), ...
            ['First input <user_mass> must be a positive scalar ' ...
            'lower than 216m/q.'])
        assert(isnumeric(timelength_mock) & ...
            2 == ndims(timelength_mock) & ...
            5 < timelength_mock & ...
            1 == size(timelength_mock,1) & ...
            1 == size(timelength_mock,2), ...
            ['Second input <timelength_mock> must be a positive ' ...
            'scalar higher than 5s.'])
        assert(ischar(varargin{1}) & ...
            2 == ndims(varargin{1}) & ...
            1 == size(varargin{1},1), ...
            ['Third input <random distribution function> must be ' ...
            'a char row vector (1xN char).'])
        assert(strcmpi(varargin{1}, 'stable') | ...
            strcmpi(varargin{1}, 'sinusoidal') | ...
            strcmpi(varargin{1}, 'constant'), ...
            ['Third input <random distribution function> must match ' ...
            'the options <stable>, <sinusoidal> or <constant>.'])
        
        start_mass = 0;
        distr_fun = varargin{1};
    otherwise
        assert(isnumeric(varargin{2}) & ...
            2 == ndims(varargin{2}) & ...
            0 < varargin{2} & ...
            1 == size(varargin{2},1) & ...
            1 == size(varargin{2},2), ...
            ['Fourth input <measured masses starting point> must ' ...
            'be a positive scalar higher than 0m/q.'])
        
        start_mass = varargin{2};
        spectrum_extension = 216; % Span of detection range of PTRMS
        end_mass_rng = start_mass + spectrum_extension;
        
        assert(isnumeric(user_mass) & ...
            2 == ndims(user_mass) & ...
            end_mass_rng > user_mass & ...
            start_mass < user_mass & ...
            1 == size(user_mass,1) & ...
            1 == size(user_mass,2), ...
            ['First input <user_mass> must be a positive scalar ' ...
            'lower than ' num2str(round(end_mass_rng)) 'm/q.'])
        assert(isnumeric(timelength_mock) & ...
            2 == ndims(timelength_mock) & ...
            5 < timelength_mock & ...
            1 == size(timelength_mock,1) & ...
            1 == size(timelength_mock,2), ...
            ['Second input <timelength_mock> must be a positive ' ...
            'scalar higher than 5s.'])
        assert(ischar(varargin{1}) & ...
            2 == ndims(varargin{1}) & ...
            1 == size(varargin{1},1), ...
            ['Third input <random distribution function> must be ' ...
            'a char row vector (1xN char).'])
        assert(strcmpi(varargin{1}, 'stable') | ...
            strcmpi(varargin{1}, 'sinusoidal') | ...
            strcmpi(varargin, 'constant'), ...
            ['Third input <random distribution function> must match ' ...
            'the options <stable>, <sinusoidal> or <constant>.'])
        
        distr_fun = varargin{1};
end


% =========================================================================
% Constants
% =========================================================================
avg_mass_step = 0.0014; % Mass resolution of the PTRMS
avg_time_step = 0.3; % Time in between ions sampling from the PTRMS
spectrum_extension = 216; % Span of detection range of PTRMS
chunk_size = 6; % Storage file chunk size for time points
cycles_pos = 27; % Position index for the beginning of the 'number of
                 % cycles' entry in the logfile of each .h5 file
shift_to_zeroes_pos = 17; % Shift from end of the 'cycles' entry to the
                          % beginning of the 'zeroes' entry
neighbourhood_mock = 0.14; % Neighbourhood of masses whose ion count still
                           % contributes to the evaluation of the central 
                           % mass value the user wants to examine. This is 
                           % not for scoping, but rather for the definition
                           % of the structure of the mockfile.


% =========================================================================
% User specified data
% =========================================================================
beg_mass_rng = start_mass + 1; % starting point of the examined mass range (starts from 
                  % 1 because it's the mass of the single proton, the
                  % smallest detectable ion
% user_mass = 154.25; % Target mass to be analysed
mass_of_choice = user_mass + 1; % We're working with the shifted mass
                                % (proton is attached, so +1 in mass)
% timelength_mock = 30; % Length of time covered by the mockfile in seconds


% =========================================================================
% Flags
% =========================================================================
no_prev_files_flag = 0; % Marks if there is already a previous mock spectra
                        % file (0 = there are previous files; 1 = there is
                        % no previous file)

                        
% =========================================================================
% Check for previous mock-files / modifies the relative flag accordingly /
% assess the file-number associated with the last mock-file available (the
% one with the highest number in the name - not the most recent one)
% =========================================================================                        
mocks = dir('PTRMSmocksequence*.h5');
n_mocks = length(mocks);

if n_mocks == 0
    n_last_mock = 0;
    no_prev_files_flag = 1;
elseif n_mocks ==1
    n_last_mock = str2double(mocks.name(end-5:end-3));
else
    ns_mocks = zeros(1,n_mocks);
    for i = 1:n_mocks
        ns_mocks(i) = str2double(mocks(i).name(end-5:end-3));
    end
    n_last_mock = max(ns_mocks);
end


% =========================================================================
% Defines name for the current mockfile to be outputted (the final file
% name has to be identified by three digits '$$$' --> max n_mockfiles is
% therefore '999')
% =========================================================================
n_curr_mock = n_last_mock + 1;

if n_curr_mock >= 100
    id_curr_mock = num2str(n_curr_mock);
elseif n_curr_mock < 100 && n_curr_mock >= 10
    id_curr_mock = ['0' num2str(n_curr_mock)];
else
    id_curr_mock = ['00' num2str(n_curr_mock)];
end

name_curr_mock = ['PTRMSmocksequence' id_curr_mock '.h5'];


% =========================================================================
% Defines the exact number of timepoints for the current mockfile being
% produced (I want each mockfile to cover approximately 30s and being ~0.3s
% the average timestep in the MS data acquisition system it means that the 
% number of timesteps is 100. As I want to be flexible and allow already
% for good adaptability towards the files the platform will be working with
% I randomly assign a number of timepoints between 96 and 102. This because
% the chunk size of the data storage system is 6. Every cycle consists of 6
% entries and if the data buffer is emptied at a point in which all 6
% points of the current cycles aren't yet acquired, the system fills the
% remaining entries with zeroes, so to have a valid array structure to
% insert inside of the .h5 file. Similarly, if the number of timepoints
% randomly generated here differs from 96 or 102, there'll be a number of
% zeroes in the final data chunk.
% =========================================================================
rng('shuffle');
low_threshold = ...
    floor((timelength_mock/avg_time_step)/chunk_size)*chunk_size;
high_threshold = low_threshold + chunk_size;
n_timepoints_curr_mock = randi([low_threshold, high_threshold]);
n_timecycles_curr_mock = ceil(n_timepoints_curr_mock/6);
remainder = rem(n_timepoints_curr_mock,chunk_size);
if remainder > 0
    n_zeroes_curr_mock = chunk_size - remainder; % From 1 to 5
else
    n_zeroes_curr_mock = 0;
end

% =========================================================================
% Creates all the .h5 structures for the current mockfile
% =========================================================================
name_time_dtst = '/TimingData/BufTimes'; % Where all the times of the
                                         % single timesteps are stored
time_size_curr_mock = [6 n_timecycles_curr_mock];
h5create(name_curr_mock,name_time_dtst,time_size_curr_mock);
%--------------------------------------------------------------------------
name_mass_dtst = '/FullSpectra/MassAxis'; % Where the vector containing all
                                          % of the measured masses is
                                          % stored
masses_length = round(spectrum_extension/avg_mass_step);
mass_size_curr_mock = [masses_length 1];
h5create(name_curr_mock,name_mass_dtst,mass_size_curr_mock);
%--------------------------------------------------------------------------
name_tofdata_dtst = '/FullSpectra/TofData'; % Where ions/s data are stored
                                            % for every time point and for
                                            % every recorded mass
tofdata_size_curr_mock = ...
    [masses_length, ...
    1, ...
    chunk_size, ...
    n_timecycles_curr_mock];
h5create(name_curr_mock,name_tofdata_dtst,tofdata_size_curr_mock);
%--------------------------------------------------------------------------
name_log_dtst = '/AcquisitionLog/Log'; % Where the info about the number of
                                       % of completed cycles and zeroes in
                                       % the last cycle is stored
% log_size_curr_mock = [1 (91 + length(num2str(n_timecycles_curr_mock)))];
log_size_curr_mock = [1 1];
h5create(name_curr_mock,name_log_dtst,log_size_curr_mock, ...
    'Datatype', 'string');


% =========================================================================
% Evaluates the starting timepoint based on the presence/absence of
% previous mockfiles
% =========================================================================
if no_prev_files_flag == 1
    start_time_curr_mock = 0;
else
    % Identifies the file to open to extract the info related to time
    if n_curr_mock >= 2 && n_curr_mock <= 10
        id_last_mock = ['00' num2str(n_last_mock)];
    elseif n_curr_mock > 10 && n_curr_mock <= 100
        id_last_mock = ['0' num2str(n_last_mock)];
    else
        id_last_mock = num2str(n_last_mock);
    end
    
    % Opens time dataset and file log
    name_last_mock = ['PTRMSmocksequence' id_last_mock '.h5'];
    data_time_last_mock = h5read(name_last_mock, name_time_dtst);
    log_last_mock = char(h5read(name_last_mock, name_log_dtst));
    
    % Extracts no. of cycles & zeroes from log
    spot = '';
    cycles = '';
    while strcmp(spot, ' ') == 0
        spot = log_last_mock(cycles_pos);
        cycles = [cycles spot];
        cycles_pos = cycles_pos + 1;
    end
    
    n_timecycles_last_mock = str2double(cycles(1:(end-1)));
    n_zeroes_last_mock = ...
        str2double(log_last_mock( ...
        cycles_pos + shift_to_zeroes_pos ...
        ));
    
    % Extracts last recorded timepoint and adds timestep to get the
    % starting point of the current mockfile
    if n_zeroes_last_mock > 0
        start_time_curr_mock = ...
            data_time_last_mock((end - n_zeroes_last_mock), end) + ...
            avg_time_step;
    else
        start_time_curr_mock = ...
            data_time_last_mock(end,end) + ...
            avg_time_step;
    end
    
end

% =========================================================================
% Switch between the functions for the creation of the mock intensity
% profile. Switch is based on the variable 'distr_fun'
% =========================================================================
switch distr_fun
%--------------------------------------------------------------------------
    case 'stable'
% =========================================================================
% Define a random distribution parameters over time of the ions count.
% This distributions defines the profile over time of the central
% mass to be analysed. It employs the probability density function 'pdf'
% based on the distribution given by 'makedist' (which in turn employs a
% 4-parameters stable distribution). Distribution parameters are randomly
% picked amongst pre-defined ranges that are found to be suitable for the
% desired profile to output.
% =========================================================================
lim_alpha = [1 2];
lim_beta = [-1 1];
lim_gamma = [10 timelength_mock];

% Defines the delta par based on the time range sampled. Delta is related 
% to the central position of the distribution, therefore it needs to 
% account for the actually measured time chunk
if no_prev_files_flag == 1
    lim_delta = [0 30];
else
    lim_delta = [(start_time_curr_mock) ...
        (start_time_curr_mock + n_timepoints_curr_mock*avg_time_step)];
end

% Randomly defines the 4 parameters within the defined ranges
lims = [lim_alpha; lim_beta; lim_gamma; lim_delta];
for i = 1:4
    if i == 1
        alpha = (lims(i,2) - lims(i,1))*rand + lims(i,1);
    elseif i == 2
        beta = (lims(i,2) - lims(i,1))*rand + lims(i,1);
    elseif i == 3
        gamma = (lims(i,2) - lims(i,1))*rand + lims(i,1);
    elseif i == 4
        delta = (lims(i,2) - lims(i,1))*rand + lims(i,1);
    end
end

% Creates central distribution
scaling_factor = 100000; % Takes output of distribution function to values
                         % similar to those given by actual ions count on
                         % real measurement of PTRMS

times = linspace(start_time_curr_mock, ...
    (start_time_curr_mock + avg_time_step*(n_timepoints_curr_mock-1)), ...
    n_timepoints_curr_mock);
distribution = ...
    makedist('Stable','alpha',alpha,'beta',beta,'gam',gamma,'delta',delta);
central_ions_distr = pdf(distribution,times);
central_ions_distr = central_ions_distr*scaling_factor;


%--------------------------------------------------------------------------
    case 'sinusoidal'
% =========================================================================
% Define a random distribution parameters over time of the ions count.
% This distributions defines the profile over time of the central
% mass to be analysed. It employs a 'sinusoidal' function, slowed down to
% guarantee a non-excessive oscillation. Distribution parameters are
% randomly picked amongst pre-defined ranges that are found to be suitable 
% for the desired profile to output.
% =========================================================================
lim_slow_factor = [2.5 4];
lim_scaling_factor = [500 5000]; % Takes output of distribution function to 
                                 % values similar to those given by actual  
                                 % ions count on real measurement of PTRMS

slow_factor = (lim_slow_factor(2) - lim_slow_factor(1))*rand + ...
    lim_slow_factor(1);
scaling_factor = (lim_scaling_factor(2) - lim_scaling_factor(1))*rand + ...
    lim_scaling_factor(1);

times = linspace(start_time_curr_mock, ...
    (start_time_curr_mock + avg_time_step*(n_timepoints_curr_mock-1)), ...
    n_timepoints_curr_mock);

central_ions_distr = scaling_factor*sin(times/slow_factor);

if no_prev_files_flag == 1
    central_ions_distr = central_ions_distr + lim_scaling_factor(2);
end


%--------------------------------------------------------------------------
    case 'constant'
% =========================================================================
% Define a random distribution parameters over time of the ions count.
% This distributions defines the profile over time of the central
% mass to be analysed. It employs a 'sinusoidal' function, slowed down to
% guarantee a non-excessive oscillation. Distribution parameters are
% randomly picked amongst pre-defined ranges that are found to be suitable 
% for the desired profile to output.
% =========================================================================
times = linspace(start_time_curr_mock, ...
    (start_time_curr_mock + avg_time_step*(n_timepoints_curr_mock-1)), ...
    n_timepoints_curr_mock);

if no_prev_files_flag == 1
    prompt = 'Output level';
    dlgtitle = 'Enter numerical value of mocked output';
    dims = [1 100];
    output = inputdlg(prompt,dlgtitle,dims);
    
    central_ions_distr = linspace( ...
        str2double(output{1}), ...
        str2double(output{1}), ...
        n_timepoints_curr_mock);
else
    central_ions_distr = zeros(1, n_timepoints_curr_mock);
end
%--------------------------------------------------------------------------
end


% =========================================================================
% Generates the orthogonal distribution related to the masses. Target mass
% isn't the only one that registers ions related to the target coumpound,
% also a small neighbourhood records ions which are fundamental to the
% accurate measurement of the concentration of the target compound.
% Therefore, it is desirable to simulate the same profile in the mock file.
% To do this, a normal distribution is generated. Combining the ions
% distribution over time with the distribution over the masses, will
% generate the surface of peaks across both masses and timepoints.
% =========================================================================
% Generates the whole list of masses that the PTRMS will analyse
data_mass_curr_mock = ...
    [beg_mass_rng:...
    avg_mass_step:...
    (beg_mass_rng + (masses_length - 1)*avg_mass_step)];
data_mass_curr_mock = data_mass_curr_mock';

% Extracts the mass range in which the mock profile will be generated
% according to the user defined mass value whose peak must be simulated
[~,ix_mass_choice] = min(abs(data_mass_curr_mock - mass_of_choice));
ixs_neighbourhood = round(neighbourhood_mock/avg_mass_step);
masses_mock = data_mass_curr_mock(...
    ix_mass_choice - ixs_neighbourhood:...
    ix_mass_choice + ixs_neighbourhood...
    );


% Generates the distribution
norm_distr = normpdf( ...
    masses_mock, ... % Created across the prev identified range of masses
    mass_of_choice, ... % Centre is on the target mass
    neighbourhood_mock/3); % Gives quite a narrow distr, okay for purpose
norm_distr = norm_distr/max(norm_distr); % Normalise to get 'height' values
                                         % of the distr to be 0:1


% =========================================================================
% Accounts for the eventual previous mockfile distribution and guarantees
% continuity across the next mockfile in the distribution of peaks. First,
% it finds the position of the last peak entry for the target mass in the
% previous mockfile, then it extracts it, and lastly it uses it to offset 
% the currently generated distribution.
% =========================================================================
if no_prev_files_flag == 0
    ix_last_entry_last_mock = chunk_size - n_zeroes_last_mock;
    read_start = [...
        ix_mass_choice, ...
        1, ...
        ix_last_entry_last_mock, ...
        n_timecycles_last_mock ...
        ];    
    read_count = [1 1 1 1];

    endvalue_tofdata_last_mock = ...
        h5read(name_last_mock,name_tofdata_dtst,read_start,read_count);

    central_ions_distr = ...
        central_ions_distr + ...
        (endvalue_tofdata_last_mock - central_ions_distr(1));
    
    if any(central_ions_distr < 0)
        central_ions_distr(central_ions_distr < 0) = 0;
    end
end


% =========================================================================
% Generates the mixed time-mass distribution of the ion peaks. Then it adds
% the zeroes to maintain proper chunk size, and lastly it structures the
% data in the same way they are stored in the original .h5 files.
% =========================================================================
cross_ions_distr = zeros(length(masses_mock),n_timepoints_curr_mock);
for i = 1:length(masses_mock)
    cross_ions_distr(i,:) = norm_distr(i)*central_ions_distr;
end

% Adds zeroes
if n_zeroes_curr_mock == 0
    cross_ions_distr_w_nulls = cross_ions_distr;
else
    cross_ions_distr_w_nulls = [ ...
        cross_ions_distr, ...
        zeros(length(masses_mock), n_zeroes_curr_mock) ...
        ];
end

% Rearrange data to be written to .h5 file
data_tofdata_curr_mock = ...
    zeros(length(masses_mock), 1, 6, n_timecycles_curr_mock);
for i = 1:length(masses_mock)
    for j = 1:n_timecycles_curr_mock
        data_tofdata_curr_mock(i,1,:,j) = ...
            cross_ions_distr_w_nulls( ...
            i, ...
            ((j-1)*chunk_size+1):((j-1)*chunk_size+chunk_size) ...
            );
    end
end


% =========================================================================
% Organises and writes data on the previously generated datasets of the
% mockfile
% =========================================================================
data_time_curr_mock = zeros(time_size_curr_mock);
if n_zeroes_curr_mock == 0
    for i = 1:n_timecycles_curr_mock
        data_time_curr_mock(:,i) = ...
            times((((i-1)*chunk_size)+1):(((i-1)*chunk_size)+chunk_size));
    end
else
    for i = 1:(n_timecycles_curr_mock - 1)
        data_time_curr_mock(:,i) = ...
            times((((i-1)*chunk_size)+1):(((i-1)*chunk_size)+chunk_size));
    end
    data_time_curr_mock(:,end) = ...
        [times(((n_timecycles_curr_mock-1)*chunk_size+1):end), ...
        zeros(1, n_zeroes_curr_mock)];
end

h5write(name_curr_mock,name_time_dtst,data_time_curr_mock);
%--------------------------------------------------------------------------
h5write(name_curr_mock,name_mass_dtst,data_mass_curr_mock);
%--------------------------------------------------------------------------
write_start = [ ...
    (ix_mass_choice - ixs_neighbourhood), ...
    1, ...
    1, ...
    1 ...
    ];
write_count = [length(masses_mock) 1 6 n_timecycles_curr_mock];

h5write( ...
    name_curr_mock, ...
    name_tofdata_dtst, ...
    data_tofdata_curr_mock, ...
    write_start, ...
    write_count ...
    );
%--------------------------------------------------------------------------
log_curr_mock = string([ ...
    'Acquisition aborted after ' ...
    num2str(n_timecycles_curr_mock) ...
    ' complete writes. ' ...
    num2str(n_zeroes_curr_mock) ...
    ' additional bufs in the incomplete last write.' ...
    ]);

h5write(name_curr_mock,name_log_dtst,log_curr_mock)

















% % =========================================================================
% % Checks if continuity of time and ions is maintained, if mockfile's
% % timespan is respected, and if log is properly written
% % =========================================================================
% if no_prev_files_flag == 0
%     bool1 = isequal(times(1),start_time_curr_mock);
%     bool2 = times(end) > (times(1) + (timelength_mock - 2)) && ...
%         times(end) < (times(1) + (timelength_mock + 2));
%     bool3 = isequal(endvalue_tofdata_last_mock,central_ions_distr(1)) |...
%         abs(endvalue_tofdata_last_mock - central_ions_distr(1)) < 1;
%     logcheck = char(log_curr_mock);
%     bool4 = strcmp(logcheck(cycles_pos:(cycles_pos + 1)), '  ');
%     
%     if bool1 == 0
%         msgbox('Time discontinuity');
%     elseif bool2 == 0
%         msgbox('Mockfile time span is out of boundaries');
%     elseif bool3 == 0
%         msgbox(['Ions profile discontinuity.' ...
%             newline 'Last value = ' num2str(endvalue_tofdata_last_mock) ...
%             newline 'First value = ' num2str(central_ions_distr(1))]);
%     elseif bool4 == 1
%         msgbox('Error in cycles number logging');
%     end
% end
% 
% 
% % =========================================================================
% % Plots the ions peaks distribution across masses and timepoints
% % =========================================================================
% figure(1)
% [timegrid, massgrid] = meshgrid(times,masses_mock);
% subplot(2,1,1);
% s = surf(timegrid, massgrid, cross_ions_distr);
% s.EdgeColor = 'none';
% 
% 
% % =========================================================================
% % Saves time and central ion profile to plot over time
% % =========================================================================
% if no_prev_files_flag == 1
%     x = times;
%     y = central_ions_distr;
% else
%     load('whole_simulation.mat')
%     x = [x times];
%     y = [y central_ions_distr];
% end
% 
% save('whole_simulation.mat','x','y')
% 
% 
% % =========================================================================
% % Plots ions peaks distribution across elapsed time since the first
% % mockfile that has been generated
% % =========================================================================
% subplot(2,1,2);
% plot(x,y)
end





