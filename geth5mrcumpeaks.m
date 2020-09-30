function [cumpeakprof, mass_rng, times] = geth5mrcumpeaks( ...
    filename, ...
    central_mass, neighbourhood)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS (or mockfile)
% 'central_mass' = Mass target to be analysed.
%                  To be considered as --> !!!(MASS + 1)!!!
% 'neighbourhood' = Defines the mass range to be included in the
%                   calculation of the peak at the specified mass target
%                   [central_mass-neighbourhood,central_mass+neighbourhood]
% 
% OUTPUTs
% 'cumpeakprof' = 1xP array containing the cumulative signal intensity
%                 profile over time (P) for the mass range defined in the
%                 inputs.
%                 Each column contains the  cumulative ions/s value
%                 calculated across all of the masses for a fixed
%                 timepoint.
% 'mass_rng' = Mrx1 column array containing the masses within the range
%              described by the central mass +- neighbourhood
% 'times' = 1xP row array containing all of the timepoints registered by
%           the current file
%                 
%
% Function to extract the cumulative signal intensity across the inputted
% mass range over time. Storage structure within .h5 file is [M,1,S,C],
% where M is the number of all of the detected masses, S is the chunk size
% (=6), and C is the number of completed acquisition cycles. Mind that this
% function extracts all of the values and rearranges them in a MrxP
% structure, where P is given by P = S*C - n_zeroes (zeroes are chunked),
% and Mr is equal to the number of mass points between the provided range
% to be analysed. Eventually, it sums up the values across each column, to
% get an output like 1xP.
% 
% DEPENDANCIES: 'geth5log.m', 'geth5mocklog.m', 'geth5masses.m',
%               'geth5time.m'  
% =========================================================================


% =========================================================================
% Initialisation and error handling
% =========================================================================
format long

assert(ischar(filename),'First input <filename> must be a char array.')

check = strcmpi(filename((end - 2):end),'.h5');
if ~check
	error('Filetype was not expected. Use .h5 file.')
end
%--------------------------------------------------------------------------
bool = isnumeric(central_mass) & isnumeric(neighbourhood);

assert(bool, ['Second and third inputs <central_mass> and ' ...
    '<neighbourhood> must be numeric values.']);
%--------------------------------------------------------------------------
% Gathers data on masses analysed to check if 'central_mass' and
% 'neighbourhood' are valid picks
masses = geth5masses(filename);

assert(~(central_mass < masses(1) | central_mass > masses(end)), ...
    ['Second input <central_mass> does not fall within the mass range ' ...
    'measured by the inputted file.'])


% =========================================================================
% Pair mass range to mass list indexes
% =========================================================================
cond = (central_mass - neighbourhood) < masses(1) | ...
    (central_mass + neighbourhood) > masses(end);

assert(~cond, 'Neighbourhood to be analysed exceeds mass range limits.')

[~,ix_min] = min(abs(masses - (central_mass - neighbourhood)));
[~,ix_max] = min(abs(masses - (central_mass + neighbourhood)));



% =========================================================================
% Log data extraction
% =========================================================================
if contains(filename, 'mock')
    [n_cycles, n_zeroes] = geth5mocklog(filename);
    chunk_size = 6; % Storage file chunk size for time points
else
    [n_cycles, n_zeroes, chunk_size, error_msg] = geth5log(filename);
    if strcmp(error_msg, 'MSmode_file')
        [n_cycles, n_zeroes, chunk_size] = geth5logMS(filename);
    end
end


% =========================================================================            
% Extracting peak data
% =========================================================================
name_tofdata_dtst = '/FullSpectra/TofData'; % Where ions/s data are stored
                                            % for every time point and for
                                            % every recorded mass
% chunk_size = 6; % Storage file chunk size for time points

start = [ix_min 1 1 1];
x_count = length(masses(ix_min:ix_max));
count = [x_count 1 chunk_size n_cycles];
stride = [1 1 1 1];
%--------------------------------------------------------------------------                                  
tofdata_raw = h5read(filename, name_tofdata_dtst, start, count, stride);

time_length = chunk_size*n_cycles;

tofdata = zeros(x_count,time_length);
for i = 1:x_count
    singlemass = zeros(1,time_length);
    for j = 1:(n_cycles)
        singlemass( ...
            (((j-1)*chunk_size)+1): ...
            (((j-1)*chunk_size)+chunk_size) ...
            ) = tofdata_raw(i,1,:,j);
    end
    tofdata(i,:) = singlemass;
end
%--------------------------------------------------------------------------
for z = 1:n_zeroes
    tofdata(:,end) = [];
end
%--------------------------------------------------------------------------
cumpeakprof = sum(tofdata);

mass_rng = masses(ix_min:ix_max);

times = geth5times(filename);

% cumpeakprof = zeros(1,length(times));
% for i = 1:length(times)
%     single = trapz(mass_rng,tofdata(:,i));
%     cumpeakprof(i) = single;
% end
end