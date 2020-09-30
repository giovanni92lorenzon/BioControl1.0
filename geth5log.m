function [n_cycles, n_zeroes, chunk_size, error_msg] = geth5log(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS in GC-mode
% 
% OUTPUTs
% 'n_cycles' = Number of completed acquisition cycles
% 'n_zeroes' = Number of zeroes in the last cycle
% 'chunk_size' = Chunk size of stored data
% 'error_msg' = Regarding unexpected n_cycles, n_zeroes, or chunk_size
%
%
% Function to extract number of completed acquisition cycles and number of
% zeroes in the final cycle FOR REGULAR PTR-MS FILES in GC MODE
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
name_log_dtst = '/AcquisitionLog/Log'; % Where the vector containing all
                                          % of the measured masses is
                                          % stored
cycles_pos = 27; % Position index for the beginning of the 'number of
                 % cycles' entry in the logfile of each .h5 file
shift_to_zeroes_pos = 17; % Shift from end of the 'cycles' entry to the
                          % beginning of the 'zeroes' entry
expected_chunk_size = 6;
%--------------------------------------------------------------------------                                   
log = h5read(filename, name_log_dtst);
log = log.logtext';

length_check = find(log(2,:) == ' ');
if length(length_check) <= 3
    error_msg = ['MSmode_file'];
    n_cycles = 0;
    n_zeroes = 0;
    chunk_size = 0;
    return
end

spot = '';
cycles = '';
while ~strcmp(spot, ' ')
    spot = log(2,cycles_pos);
    cycles = [cycles spot];
    cycles_pos = cycles_pos + 1;
end

n_cycles = str2double(cycles(1:(end-1)));


% =========================================================================
% Double check with actual time array
% =========================================================================
info = h5info(filename, '/TimingData/BufTimes');
time_dtst_size = info.Dataspace.Size;

check1 = time_dtst_size(1) ~= expected_chunk_size;
check2 = time_dtst_size(2) ~= n_cycles;

if  check1 || check2
    n_cycles = time_dtst_size(2);
    chunk_size = time_dtst_size(1);
    
    start = [1 n_cycles];
    count = [chunk_size 1];
    last_cycle_data = ...
        h5read(filename, '/TimingData/BufTimes', start, count);
    indexes = find(last_cycle_data == 0);
    n_zeroes = length(indexes);
    if check1 && ~check2
        error_msg = ['UnexpectedChunkSize(' num2str(chunk_size) ')'];
    elseif ~check1 && check2
        error_msg = ['UnexpectedCycleNumber(' num2str(n_cycles) ')'];
    else
        error_msg = ['UnexpectedChunkSize(' num2str(chunk_size) ')' ...
            '&UnexpectedCycleNumber(' num2str(n_cycles) ')'];
    end
else
    n_zeroes = str2double(log(2,(cycles_pos + shift_to_zeroes_pos)));
    chunk_size = expected_chunk_size;
    error_msg = 'None';
end