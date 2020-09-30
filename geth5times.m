function times = geth5times(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS (or mockfile)
% 
% OUTPUTs
% 'times' = 1xP row array containing all of the timepoints registered by
%           the current file
%
%
% Function to extract the timepoints measured by the PTR-MS. Storage
% structure within .h5 file is [S,C], where S is the chunk size (=6), and C
% is the number of completed acquisition cycles. Mind that the
% function extracts all of the values and rearranges them in a 1xP
% structure, where P is given by P = S*C - n_zeroes (zeroes are chunked)
% 
% DEPENDANCIES: 'geth5log.m', 'geth5mocklog.m', 'geth5logMS.m'
% =========================================================================


% =========================================================================
% Initialisation and error handling
% =========================================================================
assert(ischar(filename),'First input <filename> must be a char array.')

check = strcmpi(filename((end - 2):end),'.h5');
if ~check
	error('Filetype was not expected. Use .h5 file.')
end


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
% Time data extraction
% =========================================================================
name_time_dtst = '/TimingData/BufTimes'; % Where all the times of the
                                         % single timesteps are stored
% chunk_size = 6; % Storage file chunk size for time points
%--------------------------------------------------------------------------
time_raw = h5read(filename, name_time_dtst);

times = zeros(1,n_cycles*chunk_size);
for i = 1:(n_cycles)
    times( ...
        (((i-1)*chunk_size)+1): ...
        (((i-1)*chunk_size)+chunk_size) ...
        ) = time_raw(:,i);
end
%--------------------------------------------------------------------------
for z = 1:n_zeroes
    times(:,end) = [];
end
end