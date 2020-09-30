function tofdata = geth5allpeaks(filename)
% =========================================================================
% INPUTs
% 'filename' = Name of the .h5 file in output from the PTR-MS (or mockfile)
% 
% OUTPUTs
% 'tofdata' = MxP array containing the ions/s values recorded for all of
%             the detected masses (M) and all of the timepoints registered 
%             (P).
%             Each row contains the signal intensity profile over time for 
%             a single mass value. Each column contains the signal
%             intensity value for all of the detected masses at a single
%             timepoint.
%
%
% Function to extract all of the ions/s values over all of the masses and
% all of the timepoints. Storage structure within .h5 file is [M,1,S,C],
% where M is the number of all of the detected masses, S is the chunk size
% (=6), and C is the number of completed acquisition cycles. Mind that this
% function extracts all of the values and rearranges them in a MxP
% structure, where P is given by P = S*C - n_zeroes (zeroes are chunked).
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
% Peak intensity data extraction
% =========================================================================
name_tofdata_dtst = '/FullSpectra/TofData'; % Where ions/s data are stored
                                            % for every time point and for
                                            % every recorded mass
% chunk_size = 6; % Storage file chunk size for time points
%--------------------------------------------------------------------------                                  
tofdata_raw = h5read(filename, name_tofdata_dtst);

mass_length = size(tofdata_raw,1);
time_length = chunk_size*n_cycles;

tofdata = zeros(mass_length,time_length);
for i = 1:mass_length
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
end