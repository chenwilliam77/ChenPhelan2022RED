function out = get_ExpV_maxfnct(input, s, input_type, varargin)
% This function wraps get_ExpV into a form that is usable
% for MATLAB's various minimization routines. It also adds
% a loop over different possible thetap0 guesses to deal with
% problems satisfying the endogenous boundary conditions.
%
% Written by William Chen and Greg Phelan, Aug. 2020

% Process varargin
for j = 1:2:length(varargin)
    if strcmp(varargin{j}, 'nan_error')
        nan_error = varargin{j + 1};
    elseif strcmp(varargin{j}, 'thetap0')
        thetapgrid = varargin{j + 1};
    elseif strcmp(varargin{j}, 'min_problem')
        min_problem = varargin{j + 1};
    end
end
if ~exist('nan_error', 'var')
    nan_error = 0; % default to return the actual error
end
if ~exist('thetap0', 'var')
    thetapgrid = [s.thetap0]; % default to the thetap0 guess in s
end
if ~exist('min_problem', 'var')
    min_problem = 0; % default to maximization problem
end

% Run loop over thetapgrid
success = false;
N = length(thetapgrid);
for i = 1:N
    s.thetap0 = thetapgrid(i);
    try
        out = get_ExpV(input, s, input_type);
        success = true;
    catch err
        if i == N
            if nan_error
                out = NaN;
                return
            else
                rethrow(err);
            end
        end
    end

    if success
        break;
    else
        continue;
    end
end

% Convert to a minimization problem?
if min_problem
    out = -out;
end

end
