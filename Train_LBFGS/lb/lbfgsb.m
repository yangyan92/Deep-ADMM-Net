function [x,f,info] = lbfgsb( fcn, l, u, opts )
% x = lbfgsb( fcn, l, u )
%   uses the lbfgsb v.3.0 library (fortran files must be installed;
%       see compile_mex.m ) which is the L-BFGS-B algorithm.
%   The algorithm is similar to the L-BFGS quasi-Newton algorithm,
%   but also handles bound constraints via an active-set type iteration.
%   This version is based on the modified C code L-BFGS-B-C, and so has 
%   a slightly different calling syntax than previous versions.
%
%  The minimization problem that it solves is:
%       min_x  f(x)     subject to   l <= x <= u
%
% 'fcn' is a function handle that accepts an input, 'x',
%   and returns two outputs, 'f' (function value), and 'g' (function gradient).
%
% 'l' and 'u' are column-vectors of constraints. Set their values to Inf
%   if you want to ignore them. (You can set some values to Inf, but keep
%   others enforced).
%
% The full format of the function is:
% [x,f,info] = lbfgsb( fcn, l, u, opts )
%   where the output 'f' has the value of the function f at the final iterate
%   and 'info' is a structure with useful information
%       (self-explanatory, except for info.err. The first column of info.err
%        is the history of the function values f, and the second column
%        is the history of norm( gradient, Inf ).  )
%
%   The 'opts' structure allows you to pass further options.
%   Possible field name values:
%
%       opts.x0     The starting value (default: all zeros)
%       opts.m      Number of limited-memory vectors to use in the algorithm
%                       Try 3 <= m <= 20. (default: 5 )
%       opts.factr  Tolerance setting (see this source code for more info)
%                       (default: 1e7 ). This is later multiplied by machine epsilon
%       opts.pgtol  Another tolerance setting, relating to norm(gradient,Inf)
%                       (default: 1e-5)
%       opts.maxIts         How many iterations to allow (default: 100)
%       opts.maxTotalIts    How many iterations to allow, including linesearch iterations
%                       (default: 5000)
%       opts.printEvery     How often to display information (default: 1)
%       opts.errFcn         A function handle (or cell array of several function handles)
%                       that computes whatever you want. The output will be printed
%                       to the screen every 'printEvery' iterations. (default: [] )
%                       Results saved in columns 3 and higher of info.err variable
%
% Stephen Becker, srbecker@alumni.caltech.edu
% Feb 14, 2012
% Updated Feb 21 2015, Stephen Becker, stephen.becker@colorado.edu




narginchk(3, 4)
if nargin < 4, opts = struct([]); end

% Matlab doesn't let you use the .name convention with structures
%   if they are empty, so in that case, make the structure non-empty:
if isempty(opts), opts=struct('a',1) ; end

function out = setOpts( field, default, mn, mx )
    if ~isfield( opts, field )
        opts.(field)    = default;
    end
    out = opts.(field);
    if nargin >= 3 && ~isempty(mn) && any(out < mn), error('Value is too small'); end
    if nargin >= 4 && ~isempty(mx) && any(out > mx), error('Value is too large'); end
    opts    = rmfield( opts, field ); % so we can do a check later
end

% [f,g] = callF( x );
if iscell(fcn)
    % the user has given us separate functions to compute
    %   f (function) and g (gradient)
    callF   = @(x) fminunc_wrapper(x,fcn{1},fcn{2} );
else
    callF   = fcn;
end


n   = length(l); 
if length(u) ~= length(l), error('l and u must be same length'); end
x0  = setOpts( 'x0', zeros(n,1) );
x   = x0 + 0; % important: we want Matlab to make a copy of this. 
              %  just in case 'x' will be modified in-place
              % (Feb 2015 version of code, it should not be modified,
              %  but just-in-case, may as well leave this )
              
if size(x0,2) ~= 1, error('x0 must be a column vector'); end
if size(l,2) ~= 1, error('l must be a column vector'); end
if size(u,2) ~= 1, error('u must be a column vector'); end
if size(x,1) ~= n, error('x0 and l have mismatchig sizes'); end
if size(u,1) ~= n, error('u and l have mismatchig sizes'); end

% Number of L-BFGS memory vectors
% From the fortran driver file:
% "Values of m < 3  are not recommended, and 
%  large values of m can result in excessive computing time. 
%  The range  3 <= m <= 20 is recommended.  "
m   = setOpts( 'm', 5, 0 );


% 'nbd' is 0 if no bounds, 1 if lower bound only,
%       2 if both upper and lower bounds, and 3 if upper bound only.
% This .m file assumes l=-Inf and u=+Inf imply that there are no constraints.
% So, convert this to the fortran convention:
nbd     = isfinite(l) + isfinite(u) + 2*isinf(l).*isfinite(u);
if ispc
    nbd = int32(nbd);
else
    nbd = int64(nbd);
end


% Some scalar settings, "factr" and "pgtol"
% Their descriptions, from the fortran file:

%     factr is a DOUBLE PRECISION variable that must be set by the user.
%       It is a tolerance in the termination test for the algorithm.
%       The iteration will stop when
%
%        (f^k - f^{k+1})/max{|f^k|,|f^{k+1}|,1} <= factr*epsmch
%
%       where epsmch is the machine precision which is automatically
%       generated by the code. Typical values for factr on a computer
%       with 15 digits of accuracy in double precision are:
%       factr=1.d+12 for low accuracy;
%             1.d+7  for moderate accuracy; 
%             1.d+1  for extremely high accuracy.
%       The user can suppress this termination test by setting factr=0.
factr   = setOpts( 'factr', 1e7, 0 );

%     pgtol is a double precision variable.
%       On entry pgtol >= 0 is specified by the user.  The iteration
%         will stop when
%
%                 max{|proj g_i | i = 1, ..., n} <= pgtol
%
%         where pg_i is the ith component of the projected gradient.
%       The user can suppress this termination test by setting pgtol=0.
pgtol   = setOpts( 'pgtol', 1e-5, 0 ); % may crash if < 0

% Maximum number of outer iterations
maxIts  = setOpts( 'maxIts', 100, 1 );

% Maximum number of total iterations
%   (this includes the line search steps )
maxTotalIts     = setOpts( 'maxTotalIts', 5e3 );

% Print out information this often (and set to Inf to suppress)
printEvery  = setOpts( 'printEvery', 1 );

errFcn      = setOpts( 'errFcn', [] );

iprint  = setOpts('verbose',-1);
% <0 for no output, 0 for some, 1 for more, 99 for more, 100 for more
% I recommend you set this -1 and use the Matlab print features
% (e.g., set printEvery )

fcn_wrapper(); % initialized persistent variables
callF_wrapped = @(x,varargin) fcn_wrapper( callF, errFcn, maxIts, ...
    printEvery, x, varargin{:} );
% callF_wrapped = @(x,varargin)callF(x); % also valid, but simpler

% Call the mex file
[f,x,taskInteger,outer_count, k] = lbfgsb_wrapper( m, x, l, u, nbd, ...
    callF_wrapped, factr, pgtol, ...
    iprint, maxIts, maxTotalIts);

info.iterations     = outer_count;
info.totalIterations = k;
info.lbfgs_message1  = findTaskString( taskInteger );
errHist = fcn_wrapper();
info.err = errHist;
end % end of main function

function [f,g] = fcn_wrapper( callF, errFcn, maxIts, printEvery, x, varargin )
persistent k history
if isempty(k), k = 1; end
if nargin==0
    % reset persistent variables and return information
    if ~isempty(history) && ~isempty(k) 
        printFcn(k,history);
        f = history(1:k,:);
    end
    history = [];
    k = [];
    return;
end
if isempty( history )
    width       = 0;
    if iscell( errFcn ), width = length(errFcn);
    elseif ~isempty(errFcn), width = 1; end
    width       = width + 2; % include fcn and norm(grad) as well
    history     = zeros( maxIts, width );
end

% Find function value and gradient:
[f,g] = callF(x);

if nargin > 5
    outerIter = varargin{1}+1;
    
    history(outerIter,1)    = f;
    history(outerIter,2)    = norm(g,Inf); % g is not projected
    if isa( errFcn, 'function_handle' )
        history(outerIter,3) = errFcn(x);
    elseif iscell( errFcn )
        for j = 1:length(errFcn)
            history(outer_count,j+2) = errFcn{j}(x);
        end
    end
    
    if outerIter > k
        % Display info from *previous* input
        % Since this may be called several times before outerIter
        % is actually updated
%         fprintf('At iterate %5d, f(x)= %.2e, ||grad||_infty = %.2e [MATLAB]\n',...
%             k,history(k,1),history(k,2) );
        if ~isinf(printEvery) && ~mod(k,printEvery)
            printFcn(k,history);
            net = weiTOnet(x);
            err = history(k, 1);
            save(strcat('./Train_output/net/net-', saveName(k, 3), '.mat'), 'net');
            save(strcat('./Train_output/error/error-', saveName(k, 3), '.mat'), 'err');
        end
        k = outerIter;
    end

    
end

end


function printFcn(k,history)
fprintf('Iter %5d, loss = %f, ||grad||_infty = %f', ...
    k, history(k,1),  history(k,2) );
for col = 3:size(history,2)
    fprintf(', %.2e', history(k,col) );
end
fprintf('\n');
end
    


function [f,g] = fminunc_wrapper(x,F,G)
% [f,g] = fminunc_wrapper( x, F, G )
%   for use with Matlab's "fminunc"
f = F(x);
if nargin > 2 && nargout > 1
    g = G(x);
end

end




function str = findTaskString( taskInteger )
% See the #define statements in lbfgsb.h
switch taskInteger
case 209
    str = 'ERROR: N .LE. 0';
case 210
    str = 'ERROR: M .LE. 0';
case 211
    str = 'ERROR: FACTR .LT. 0';
case 3
	str = 'ABNORMAL_TERMINATION_IN_LNSRCH.';
case 4
	str = 'RESTART_FROM_LNSRCH.';
case 21
	str = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL.';
case 22
	str = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH.';
case 31
	str = 'STOP: CPU EXCEEDING THE TIME LIMIT.';
case 32
	str = 'STOP: TOTAL NO. of f AND g EVALUATIONS EXCEEDS LIM.';
case 33
	str = 'STOP: THE PROJECTED GRADIENT IS SUFFICIENTLY SMALL.';
case 101
	str = 'WARNING: ROUNDING ERRORS PREVENT PROGRESS';
case 102
	str = 'WARNING: XTOL TEST SATISIED';
case 103
	str = 'WARNING: STP = STPMAX';
case 104
	str = 'WARNING: STP = STPMIN';
case 201
	str = 'ERROR: STP .LT. STPMIN';
case 202
	str = 'ERROR: STP .GT. STPMAX';
case 203
	str = 'ERROR: INITIAL G .GE. ZERO ';
case 204
	str = 'ERROR: FTOL .LT. ZERO';
case 205
	str = 'ERROR: GTOL .LT. ZERO';
case 206
	str = 'ERROR: XTOL .LT. ZERO';
case 207
	str = 'ERROR: STPMIN .LT. ZERO';
case 208
	str = 'ERROR: STPMAX .LT. STPMIN';
case 212
	str = 'ERROR: INVALID NBD';
case 213
	str = 'ERROR: NO FEASIBLE SOLUTION';
    otherwise
        str = 'UNRECOGNIZED EXIT FLAG';
end
end
