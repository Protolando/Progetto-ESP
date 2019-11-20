function [y] = DifferentialEvolution(f, ndim, sampleSize, niter, a, b)
%Initialization
x=linspace(a,b,sampleSize);
x=repmat(x,ndim,1);

for i = 1:niter
    %calculate F
    F = 2-(i/niter)*2;

    %initialize the next generation vector
    nextGen = zeros(ndim, sampleSize);
    
    for j = 1:sampleSize
        n = sampleSize;
        l = 3;
        %Choose 3 random members of the population
        Xs = randperm(n,l);
        
        %calculate Xc
        Xc = transpose(x(:, Xs(1)) + F*(x(:, Xs(2))-x(:, Xs(3))));
        %calculate X'i
        nextGen(:, j) = crossover(x(:, j), Xc, ndim);
        
        %Make sure it's within the problem bounds
        for k = 1:ndim
            if(nextGen(k, j) < a)
                nextGen(k, j) = a;
            elseif(nextGen(k, j) > b)
                nextGen(k, j) = b;
            end
        end
    end
    
    %Only keep those bigger then their old generation counterpart
    for j = 1:sampleSize
        y1 = f(x(:, j));
        y2 = f(nextGen(:, j));
        if(y1 > y2)
            x(:, j) = x(:, j);
        else
            x(:, j) = nextGen(:, j);
        end
    end
end

%return the best value.
[~, i] = max(f(x));
x = x(:, i);
y = f(x(:));
end

function x = crossover(x1, x2, ndim)
    %Crossover probability
    CR = 0.5;
    
    x = x1;
    %Find all indexes for which a random number is smaller then the
    %crossover probability
    idxs = find(rand(1, ndim) < CR);
    %swap said indexes
    if(~isempty(idxs))
        x(idxs) = x2(idxs);
    end
end