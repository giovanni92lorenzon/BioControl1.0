function [n_cycles, n_zeroes] = geth5mocklog(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 mockfile to simulate live PTR-MS acquisition
% 
% OUTPUTs
% 'n_cycles' = Number of completed acquisition cycles
% 'n_zeroes' = Number of zeroes in the last cycle
%
%
% Function to extract number of completed acquisition cycles and number of
% zeroes in the final cycle FOR MOCK PTR-MS FILES
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
%--------------------------------------------------------------------------                                   
log = char(h5read(filename, name_log_dtst));

spot = '';
cycles = '';
while ~strcmp(spot, ' ')
    spot = log(cycles_pos);
    cycles = [cycles spot];
    cycles_pos = cycles_pos + 1;
end

n_cycles = str2double(cycles(1:(end-1)));
n_zeroes = str2double(log(cycles_pos + shift_to_zeroes_pos));
end