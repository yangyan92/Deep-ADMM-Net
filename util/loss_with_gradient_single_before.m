function [ loss, real_x ] = loss_with_gradient_single_before( data, net )
train = data.train;
label = data.label;
LL = double(-1:0.02:1);
N = numel(net.layers);

res = struct(...
    'x',cell(1,N+1),...
    'dzdx',cell(1,N+1),...
    'dzdw',cell(1,N+1));
res(1).x = train;

%% forword propagation 
for n = 1:N
    l = net.layers{n};
    switch l.type
        case 'X_org'
            res(n+1).x = xorg (res(n).x , l.weights{1} , l.weights{2});
        case 'Convo'
            res(n+1).x = convo(res(n).x , l.weights{1});
        case 'Non_linorg'
          res(n+1).x = zorg( LL ,res(n).x , l.weights{1});
        case 'Multi_org'
            res(n+1).x = betaorg(res(n-1).x , res(n).x ,l.weights{1} );
        case 'X_mid'
            res(n+1).x = xmid( res(n-1).x , res(n).x, res(1).x , l.weights{1} , l.weights{2});
        case 'Non_linmid'
           res(n+1).x = zmid(LL, res(n).x , res(n-2).x , l.weights{1} ); 
        case 'Multi_mid'
            res(n+1).x = betamid(res(n-3).x , res(n-1).x , res(n).x ,l.weights{1});
        case 'Multi_final'
            res(n+1).x = betafinal(res(n-3).x , res(n-1).x , res(n).x, l.weights{1});
        case 'loss'
              res(n+1).x = rnnloss(res(n).x, label); 
        otherwise
            error('No such layers type.');
      end;
    
end;
    loss = res(end).x;
    loss = double(loss);
    real_x = res(end-1).x;      
end

    
            










