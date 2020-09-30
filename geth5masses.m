function masses = geth5masses(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS (or mockfile)
% 
% OUTPUTs
% 'masses' = Mx1 column array containing all of the masses detected by the
%            PTR-MS (or mock masses)
%
%
% Function to extract the whole array of the masses analysed by the PTR-MS
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
% Mass data extraction
% =========================================================================
name_mass_dtst = '/FullSpectra/MassAxis'; % Where the vector containing all
                                          % of the measured masses is
                                          % stored
%--------------------------------------------------------------------------
masses = h5read(filename, name_mass_dtst);
end

