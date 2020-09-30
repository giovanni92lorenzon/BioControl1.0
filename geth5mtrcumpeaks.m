function [cumpeakprofrng, mass_rng, times_rng] = geth5mtrcumpeaks( ...
    filename, ...
    central_mass, neighbourhood, ...
    t_start, t_end)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS (or mockfile)
% 'central_mass' = Mass target to be analysed
% 'neighbourhood' = Defines the mass range to be included in the
%                   calculation of the peak at the specified mass target
%                   [central_mass-neighbourhood,central_mass+neighbourhood]
% 't_start' = Time (in seconds) to the define the starting point of the
%             output array
% 't_end' = Time (in seconds) to the define the end point of the output
%           array
% 
% OUTPUTs
% 'cumpeakprofrng' = 1xPr array containing the cumulative signal intensity
%                    profile over a specific range of time (Pr) within that
%                    recorded by the file in analysis for the mass range
%                    defined in the inputs.
%                    Each column contains the  cumulative ions/s value
%                    calculated across all of the masses for a fixed
%                    timepoint.
% 'mass_rng' = Mrx1 column array containing the masses within the range
%              described by the central mass +- neighbourhood
% 'times_rng' = 1xPr row array containing the timepoints within the range
%               given by the last two inputs
%                 
%
% Function to extract the cumulative signal intensity across the inputted
% mass range over a specific tange of time. Storage structure within .h5
% file is [M,1,S,C], where M is the number of all of the detected masses,
% S is the chunk size (=6), and C is the number of completed acquisition
% cycles. Mind that this function extracts all of the values and rearranges
% them in a MrxPr structure, where Pr is the selected range of times 
% belonging to the original set given by P = S*C - n_zeroes (zeroes are 
% chunked), and Mr is equal to the number of mass points between the 
% provided range to be analysed. Eventually, it sums up the values across 
% each column, to get an output like 1xPr.
% 
% DEPENDANCIES: 'geth5masses.m', 'geth5time.m', 'geth5log.m'          
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
%--------------------------------------------------------------------------
% Gathers data on times analysed to check if 't_start' and 't_end'
% are valid picks
times = geth5times(filename);

assert(~(t_start < times(1) | t_end > times(end)), ...
    'Inputted time limits exceed those of the file in exam.')


% =========================================================================
% Pair mass range to mass list indexes
% =========================================================================
cond = (central_mass - neighbourhood) < masses(1) | ...
    (central_mass + neighbourhood) > masses(end);

assert(~cond, 'Neighbourhood to be analysed exceeds mass range limits')

[~,ix_min] = min(abs(masses - (central_mass - neighbourhood)));
[~,ix_max] = min(abs(masses - (central_mass + neighbourhood)));


% =========================================================================
% Extracts chunk size to be employed
% =========================================================================
if contains(filename, 'mock')
    chunk_size = 6;
else
    [~,~,chunk_size, error_msg] = geth5log(filename);
    if strcmp(error_msg, 'MSmode_file')
        [~, ~, chunk_size] = geth5logMS(filename);
    end
end


% =========================================================================
% Pair time range to times list indexes
% =========================================================================
[~,ix_st] = min(abs(times - t_start),[],'all','linear');
[~,ix_end] = min(abs(times - t_end),[],'all','linear');

cyc_st = ceil(ix_st/chunk_size);
rem_st = rem(ix_st,chunk_size);
cyc_end = ceil(ix_end/chunk_size);
rem_end = rem(ix_end,chunk_size);
rng = cyc_end - cyc_st;

times_rng = times(ix_st:ix_end);


% =========================================================================            
% Extracting peak data
% =========================================================================
name_tofdata_dtst = '/FullSpectra/TofData'; % Where ions/s data are stored
                                            % for every time point and for
                                            % every recorded mass

start = [ix_min 1 1 cyc_st];
x_count = length(masses(ix_min:ix_max));
z_count = rng + 1;
count = [x_count 1 chunk_size z_count];
stride = [1 1 1 1];
%--------------------------------------------------------------------------                                  
tofdata_raw = h5read(filename, name_tofdata_dtst, start, count, stride);

time_length = chunk_size*z_count;

tofdata = zeros(x_count,time_length);
for i = 1:x_count
    singlemass = zeros(1,time_length);
    for j = 1:(z_count)
        singlemass( ...
            (((j-1)*chunk_size)+1): ...
            (((j-1)*chunk_size)+chunk_size) ...
            ) = tofdata_raw(i,1,:,j);
    end
    tofdata(i,:) = singlemass;
end
%--------------------------------------------------------------------------
if rem_st == 0
    if rem_end == 0
        tofdata = tofdata(1:end,chunk_size:end);
    else
        tofdata = tofdata(1:end,chunk_size:(end-(chunk_size-rem_end)));
    end
else
    if rem_end == 0
        tofdata = tofdata(1:end,rem_st:end);
    else
        tofdata = tofdata(1:end,rem_st:(end-(chunk_size-rem_end)));
    end
end
%--------------------------------------------------------------------------
cumpeakprofrng = sum(tofdata);

mass_rng = masses(ix_min:ix_max);
end