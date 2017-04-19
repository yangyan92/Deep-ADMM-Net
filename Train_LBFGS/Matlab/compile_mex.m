%{
Run this file to compile the mex files if you need to
(or if the binaries work for your computer, then you can skip this)
This file also contains two simple test problems

Stephen Becker, Feb 14, 2012  srbecker@alumni.caltech.edu
Updated         May  3, 2012  adding the f2c version for Windows userss
Updated         Jun 24, 2012  adding a lot of trouble-shooting support mainly for 64-bit linux
MAJOR Update    Feb 21, 2015  stephen.becker@colorado.edu
    completely redid mex file and lbfgs library
    (converted LBFGSB fortran code to C). Many changes. Most of the Matlab
    interface remains the same, but mex interface has changed.
    Should not affect end-user. Compilation is SO MUCH EASIER than it
    was with Fortran
    Most of the C code is the same as the Fortran, but a few fcn signatures
    different, and chanted "task" from str to int, and modified print functions
    No longer able to print to file 

See also lbfgsb.m and lbfgsb_wrapper.c
%}



%% -- Compile --

% run mex -setup if you haven't already...

SRC_DIR = fullfile('..','src');
% Note: make sure to run be in the directory containing this file
% when you run this, otherwise above path is wrong
INCLUDES = SRC_DIR;
SRC = {'lbfgsb.c','linesearch.c','subalgorithms.c','print.c',...
    'linpack.c','miniCBLAS.c','timer.c'};
for i = 1:length(SRC)
    SRC{i} = fullfile(SRC_DIR,SRC{i});
end

if ispc
    % do not specify -lm flag
    mex('lbfgsb_wrapper.c','-largeArrayDims','-UDEBUG',...
        ['-I',SRC_DIR], SRC{:} );
else
    mex('lbfgsb_wrapper.c','-largeArrayDims','-lm','-UDEBUG',...
        ['-I',SRC_DIR], SRC{:} );
end

%% test the new function
disp('=== lbfgsb "driver1" test problem (Rosenbrock, 25 dimensions) === ');
% Here's the test problem included with lbfgsb called 'driver1'
% (It's a version of the Rosenbrock test function)

n   = 25;

l   = ones(n,1); u = l;
odd = 1:2:n;
even= 2:2:n;
l(odd) = 1.0;
u(odd) = 1.0e2;
l(even)= -1.0e2;
u(even)=  1.0e2;

opts    = struct( 'x0', 3*ones(n,1) );
opts.printEvery     = 2; % controls how often we print output from .m file wrapper
opts.m  = 5;
%opts.maxIts = 10;
opts.errFcn = @(x) norm(x-1); % for now just an arbitrary fcn
% opts.verbose = -1; % default is -1, i.e., no output from mex

[x,f,info] = lbfgsb( @driver1, l, u, opts );

% The true objective value is 0.
if abs(f) < 1e-8
    disp('Success!');
    semilogy( abs(info.err(:,1)-f),'o-' ); 
    xlabel('iteration'); 
    ylabel('error in objective function');
else
    disp('Something didn''t work right :-(  ');
end

% the structure info.err contains the objective function (1st column)
%   and norm(gradient,Inf) (2nd column)


%% another test function, the 2D Rosenbrock function
disp('=== Rosenbrock test function, 2D === ');
n = 2;

fxy = @(x,y) 100*( y-x.^2).^2  +  (1-x ).^2 ;
f   = @(x)   fxy( x(1,:), x(2,:) );
gxy = @(x,y) [100*(4*x.^3-4*x.*y)+2*x-2; 100*(2*y-2*x.^2)];
g   = @(x)   gxy( x(1,:), x(2,:) );

% There are no constraints
l   = -inf(n,1);
u   = inf(n,1);

opts    = struct( 'x0', [-1.9;2] );
opts.printEvery     = 1;
opts.m  = 5;

% Here's an example of using an error function. For Rosenbrock,
%   we know the true solution, so we can measure the error at every
%   iteration:
trueSoln = [1;1];
% "errFcn" will be printed to the screen
opts.errFcn     = @(x) norm(x-trueSoln)/max(norm(trueSoln),1);
% "outputFcn" will save values in the "info" output
opts.outputFcn  = opts.errFcn;

% Ask for very high accuracy
opts.pgtol      = 1e-10;
opts.factr      = 1e3;

% The {f,g} is another way to call it
[x,f,info] = lbfgsb( {f,g} , l, u, opts );

if abs(f) < 1e-8
    disp('Success!');
% since we included opts.outputFcn, the info.err now has 3 columns.
%   The first 2 columns are the same as before; the 3rd column
%   is the output of our outputFcn
semilogy( info.err(:,3)-f,'o-' ); xlabel('iteration'); ylabel('relative error in iterate function');
else
    disp('Something didn''t work right :-(  ');
end

