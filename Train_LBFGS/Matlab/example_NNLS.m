%% Non-negative least-squares (NNLS) using L-BFBS-B 
%
% Non-negative least-squares solves the following problem:
%   $$ \min_x \|Ax \textrm{--} b\|^2_2 \quad\textrm{such that}\quad x \ge 0 $$
%
% The matrix 'A' may have more columns than rows (the 'underdetermined' case),
%   or more rows than columns (the 'overdetermined' case), or the same
%   number of rows and columns. Some solvers, such as the PQN
%   method described in "Tackling Box-Constrained Optimization via a new
%   Projected Quasi-Newton Approach" by Dongmin Kim, Suvrit Sra, and
%   Inderjit Dhillon (http://www.cs.utexas.edu/users/inderjit/public_papers/pqnj_sisc10.pdf),
%   only work for the overdetermined case.
%
% To quote from that paper, 
% "Not surprisingly, some constrained
% optimization methods have also been applied to solve NNLS. It is interesting
% to note that for large scale problems these specialized algorithms are outperformed
% by modern methods such as TRON, LBFGS-B, or the methods of this paper. Curiously
% this fact has not yet been widely adopted by the wider research community
% (footnote: This could be because Matlab continues to ship the antiquated 
%   lsqnonneg function, which is an implementation of the original NNLS algorithm of
%   Lawson and Hanson 1974 )."
%
% The Kim/Sra/Dhillon paper compares the following algorithms:
%
% Fast NNLS, by Rasmus Bro. Available at:
%   http://www.mathworks.com/matlabcentral/fileexchange/3388-nnls-and-constrained-regression
%
% mtron, mex wrapper by Christoph Ortner, available at:
%   http://www.mathworks.com/matlabcentral/fileexchange/14848-mtron
%   Based on the fortran tron algorithm by Chih-Jen Lin and Jorge More, 
%   "Newton's method for large bound-constrained optimization problems",
%   SIAM Journal on Optimization, 9(4), pp. 1100-1127, 1999.
%        http://www-unix.mcs.anl.gov/~more/tron/ 
%
% L-BFGS-B.
%   R. Byrd, P. Lu, J. Nocedal, and C. Zhu, "A Limited Memory Algorithm
% for Bound Constrained Optimization", SIAM Journal on Scientific Computing, 16
% (1995), pp. 1190--1208.
%
%
% For a published version of this demo, see
%   http://www.mathworks.com/examples/matlab/4139-non-negative-least-squares-nnls-using-l-bfbs-b
%% This demo
% Here, we use the mex wrapper for L-BFGS-B v3.0, which is a significantly
%   improved version of L-BFGS-B from v2.1. We show how to use
%   the software and the fminunc_wrapper helper file.
%
% It also compares to some NNLS implementations availabe on the matworks
%   file exchange. In addition to Fast NNLS (FNNLS), mtron, and LBFGS,
%   we compare with the following algorithms, all written by Uriel Roque
%   and based on: Portugal, Judice and Vicente, 
% "A comparison of block pivoting and interior point algorithms for 
% linear least squares problems with nonnegative variables",
%  Mathematics of Computation, 63(1994), pp. 625-643
%
% activeset.m   This is pretty fast for medium-scale and smaller problems
%       http://www.mathworks.com/matlabcentral/fileexchange/10908-active-set-algorithm
%
% blocknnls.m   Similar to activeset.m in performance
%       http://www.mathworks.com/matlabcentral/fileexchange/8157-nnls/content/blocknnls.m
%
% newton.m  Very slow for large problems
%       http://www.mathworks.com/matlabcentral/fileexchange/10953-newton-s-algorithm-for-nnls
% 
% pcnnls.m (predictor-corrector method) Very slow for large problems
%       http://www.mathworks.com/matlabcentral/fileexchange/8150-predictor-corrector-algorithm
%
% The most interesting tests use large matrices. For small matrices, tests
% are pointless, because any of the methods are suitable.
%% Setup a problem

% The best codes handle N = 20,000 as long as the matrix is very sparse.
% N   = 3000; M = 4000; % Large scale. Things start to get interesting
N   = 1000; M = 1500;     % at this size, some algo take a long time!
% N   = 100; M = 150;     % at this size, all algorithms take < 14 seconds
A   = randn(M,N);
b   = randn(M,1);

fcn     = @(x) norm( A*x - b)^2;
% here are two equivalent ways to make the gradient. grad2 is sometimes faster
grad1    = @(x) 2*A'*(A*x-b);
AtA     = A'*A; Ab = A'*b;
grad2    = @(x) 2*( AtA*x - Ab );

grad    = grad2;


x = [];
time = [];

%% Solve NNLS with L-BFGS-B

l  = zeros(N,1);    % lower bound
u  = inf(N,1);      % there is no upper bound
tstart=tic;
fun     = @(x)fminunc_wrapper( x, fcn, grad); 
% Request very high accuracy for this test:
opts    = struct( 'factr', 1e4, 'pgtol', 1e-8, 'm', 10);
opts.printEvery     = 5;
if N > 10000
    opts.m  = 50;
end
% Run the algorithm:
[xk, ~, info] = lbfgsb(fun, l, u, opts );
t=toc(tstart)
% Record results
x.lbfgsb    = xk;
time.lbfgsb = t;


%% Solve with TRON, via MTRON interface

% Only run this if you have mtron installed and it is in the path
if exist( 'itron.m', 'file' )

    x0   = zeros(N,1);
    xl   = zeros(N,1);
    xu   = +1e300*ones(N,1);
    fmin = -1e300;
    H       = sparse(AtA/2); % will crash if not a sparse matrix
    tstart=tic;
    hess    = @(x) H;
    fun     = @(x)fminunc_wrapper( x, fcn, grad, hess );
    [xk, fval, exitflag, output] = itron(fun, x0, xl, xu, fmin );
    t=toc(tstart)
    x.tron    = xk;
    time.tron = t;
end

%% Active set. Fast on medium problems
if exist( 'activeset.m', 'file' )
    tstart=tic;
    [xk,y]  = activeset(A,b);
    t=toc(tstart)
    x.activeset    = xk;
    time.activeset = t;
end
%% Block pivoting. Fast on medium problems
if exist( 'blocknnls.m', 'file' )
    tstart=tic;
    [xk]  = blocknnls(A,b, 'fixed');
    t=toc(tstart)
    x.blockPivot    = xk;
    time.blockPivot = t;
end
%% Newton.  Slow!
if exist( 'newton.m', 'file' ) && N < 500
    tstart=tic;
    [xk,y]  = newton(A,b, ones(N,1), 100); % can't have 0 starting vector
    t=toc(tstart)
    x.newton    = xk;
    time.newton = t;
else
    fprintf('Skipping Newton method because we can''t find it, or it is too slow\n');
end
%% Predictor-Corrector. Can be very slow
if exist( 'pcnnls.m', 'file' ) && N < 500
    tstart=tic;
    [xk,y,nits]  = pcnnls(A,b,ones(N,1), 3000);
    t=toc(tstart)
    x.predCorr    = xk;
    time.predCorr = t;
else
    fprintf('Skipping predCorr method because we can''t find it, or it is too slow\n');
end
%% Run Matlab's default (Lawson and Hanson) Very slow on large problems
tstart=tic;
xk = lsqnonneg(A,b);
t=toc(tstart)
x.lsqnonneg   = xk;
time.lsqnonneg = t;
%% Fast NNLS, modification of Lawson and Hanson. Much better for large problems
if exist( 'fnnls.m', 'file' )
    tstart=tic;
    [xk]  = fnnls(A'*A,A'*b);
    t=toc(tstart)
    x.fnnls    = xk;
    time.fnnls = t;
end
%% PQN-LBFGS and PQN-BB algorithms of Kim/Sra/Dhillon. Very fast.
if exist( 'solnls.m', 'file' )
    opt = solopt;
    opt.maxtime     = 2000;
    opt.verbose     = 0;
    tstart=tic;
    % run their 'BB' variant
    opt.algo = 'BB';
    out     = solnls( A, b, zeros(N,1), opt );
    t=toc(tstart)
    x.PQN_BB   = out.x;
    time.PQN_BB = t;

    % and run their 'PLB' variant (their 'PQN' variant is much slower)
    %   which uses L-BFGS (not to be confused with L-BFGS-B)
    opt.algo = 'PLB';

    tstart=tic;
    out     = solnls( A, b, zeros(N,1), opt );
    t=toc(tstart)

    x.PQN   = out.x;
    time.PQN = t;
    
end
%% Results
% Find the best answer, and use that as the reference.
fMin = Inf;
for f=fieldnames(x)',
    if fcn(x.(f{1})) < fMin,
        fMin = fcn(x.(f{1}));
        best = f{1};
    end
end
xReference = x.(best);
errFcn      = @(x) norm(x-xReference)/norm(xReference);

% Print out info. Verify that the solution is indeed non-negative (hence the
%   min(x) information), and the objective function, and the error
%   against the reference solution. Also display the time.
fprintf('== Size of problem is %d x %d == \n', M, N );
for f=fieldnames(x)',
    fprintf('%10s:  obj is %7.2f, min(x) is %7.1d, err is %.2e, time is %6.3f s\n', ...
        f{1}, fcn(x.(f{1})), min(x.(f{1})), errFcn(x.(f{1})), time.(f{1}) );
end
