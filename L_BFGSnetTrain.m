%% This is a traning code for ADMM-Net by L-BFGS optimizing.
%% If you use this code, please cite our paper:
%% [1] Yan Yang, Jian Sun, Huibin Li, Zongben Xu. Deep ADMM-Net for Compressive Sensing MRI, NIPS(2016).
%% Copyright (c) 2017 Yan Yang
%% All rights reserved.

  clear all ;
  clc;

%% Network initialization
net = InitNet ( );

%% Initial loss
wei0 = netTOwei(net);
l0 = loss_with_gradient_total(wei0)

%% L-BFGS optimiztion
fun = @loss_with_gradient_total;
%parameters in the L-BFGS algorithm
low = -inf*ones(length(wei0),1);
upp = inf*ones(length(wei0),1);
opts.x0 = double(gather(wei0));
opts.m = 5;
opts.maxIts = 7.2e4;
opts.maxTotalIts = 7.2e4;
opts.printEvery = 1;
[wei1, l1, info] = lbfgsb(fun, low, upp, opts);
wei1=single(wei1);
net1 = weiTOnet(wei1);
fprintf('Before training, error is %f; after training, error is %f.\n', l0, l1);
