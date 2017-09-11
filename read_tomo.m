function tomo = read_tomo(filename)
%
% USEAGE: TOMO = read_tomo(FILENAME)
%
% INPUT ARGUMENTS:
%
% FILENAME:
%  The filename of the tomograph-parameter file.
%
% OUTPUTS:
%
% TOMO
%  A 1xm struct array containg the tomograph parameters. Fields include:
%       FIELD       DESCRIPTION
%       filename    The ring-parameter filename.
%
%       shift       The axial position of the ring.
%
%       rot         The transaxial rotation, a counterclockwise rotation
%                   between 0 and 360 degrees.

if (exist(filename, 'file') == 2)
    tomo_params = fileread(filename);
else
    error('The tomograph-parameter file "%s" cannot be found.', filename);
end

% filename
fn_expr = 'STR\s+blocktomo_ring_parameter_file\s*=\s*"([a-zA-Z_-.]+)"';
fn = regexp(tomo_params, fn_expr, 'tokens');

% shift
shift_str = 'REAL\s+blocktomo_ring_axial_shift\s*=\s*';
shift = get_val(shift_str, tomo_params);

% rot
rot_str = 'REAL\s+blocktomo_ring_transaxial_rotation\s*=\s*';
rot = get_val(rot_str, tomo_params);

% Display an error if the number of each parameter read isn't the same
num_rings = numel(fn);
if ( num_rings ~= [numel(shift), numel(rot)])
    error('The parameter file is not formatted correctly.');
end

for i = 1:num_rings
    tomo(i) = struct(   'filename', fn{i},      ...
                        'shift',    shift(i),   ...
                        'rot',      rot(i) );
end

end

function val = get_val(pattern, params)
% Parses params for the given pattern which is a line which ends with a value.
%
% USAGE: VAL = get_val(PATTERN, PARAMS)
%
% INPUT ARGUMENTS:
%
% PATTERN
%  A character array which contains the beginning of the line which ends
%  with the value which is desired.
%
% PARAMS
%  A character vector of the block-parameter file.
%
% OUTPUT: VAL
%  A 1xm vector with the values embedded in the pattern, where m is the number
%  of lines that match in the pattern in params.

expr = [pattern, '(-?[0-9]+\.?[0-9]*)'];
out = regexp(params, expr, 'tokens');
if (~isempty(out))
    for m = 1:numel(out)
        val(m) = str2double(out{m});
    end
else
    error([ 'The parameter file does not contain the value ', ...
            'in the proper format.']);
end

end
