function [n_cycles, n_zeroes, chunk_size] = geth5logMS(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS in MS-mode
% 
% OUTPUTs
% 'n_cycles' = Number of completed acquisition cycles
% 'n_zeroes' = Number of zeroes in the last cycle
% 'chunk_size' = Chunk size of stored data
%
%
% Function to extract number of completed acquisition cycles and number of
% zeroes in the final cycle FOR REGULAR PTR-MS FILES in MS MODE
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
info = h5info(filename, '/TimingData/BufTimes');
time_dtst_size = info.Dataspace.Size;

n_cycles = time_dtst_size(2);
chunk_size = time_dtst_size(1);

start = [1 n_cycles];
count = [chunk_size 1];
last_cycle_data = ...
    h5read(filename, '/TimingData/BufTimes', start, count);
indexes = find(last_cycle_data == 0);
n_zeroes = length(indexes);

end