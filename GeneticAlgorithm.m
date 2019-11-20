function [y] = GeneticAlgorithm(f, ndim, populationSize, niter, a, b)

%Initialization
digitsBefComma = 3;
digitsAftComma = 10;
nbits = 1 + digitsBefComma + digitsAftComma; %sign + integer digits + decimal digits

%generate the population matrix with linspace
x = zeros(ndim, nbits, populationSize);
for k = 1:populationSize
    v = rand(1, ndim).*(b-a)+a;
    x(:, :, k) = float2bin(v(:));
end

% %initialize the results
y = zeros(1, populationSize);
for i = 1:populationSize
    y(i) = f(bin2float(x(:, :, i)));
end

%sort according to ranking
[~, order] = sort(y, 'descend');
x = x(:, :, order);
y = y(order);

%calculate probabilities for this iteration
beta = 2;%*i/niter;
probabilities = (1/(populationSize))*(beta-2*(beta-1)*(linspace(1,populationSize, populationSize)-1)/(populationSize-1));
%calculate PMF(Probability Mass Function)
probabilities = cumsum(probabilities);

for i = 1:niter
    %Initialize the new generation array and result
    nextGenX = zeros(ndim, nbits, populationSize);
    nextGenY = zeros(1, populationSize);
    
    for j = 1:2:populationSize
        for h = 1:ndim
            %choose randomly 2 vectors according to the probability
            v = [find(probabilities>rand(1),1) find(probabilities>rand(1),1)];
            %Crossover for every dimention
            [nextGenX(h,:,j),nextGenX(h,:,j+1)] = crossover(x(h, :, v(1)), x(h, :, v(2)));
        end
        
        %Mutate both vectors
        nextGenX(:,:,j) = mutate(nextGenX(:,:,j));
        %Ensure the value is within the problem bounds
        val = bin2float(nextGenX(:,:,j));
        if(val < a)
            nextGenX(:,:,j) = float2bin(a);
        elseif(val > b)
            nextGenX(:,:,j) = float2bin(b);
        end

        nextGenX(:,:,j+1) = mutate(nextGenX(:,:,j+1));
        val = bin2float(nextGenX(:,:,j+1));
        if(val < a)
            nextGenX(:,:,j+1) = float2bin(a);
        elseif(val > b)
            nextGenX(:,:,j+1) = float2bin(b);
        end
        
        nextGenY(j) = f(bin2float(nextGenX(:,:,j)));
        nextGenY(j+1) = f(bin2float(nextGenX(:,:,j+1)));
    end
    
    maxnew = max(nextGenY);
    %initialize an empty temporary matrix
    newy = zeros(1, populationSize);
    newx = zeros(ndim, nbits, populationSize);

    [~, order] = sort(nextGenY, 'descend');
    nextGenX = nextGenX(:, :, order);
    nextGenY = nextGenY(order);
    
    %Elitism
    val = 0;
    for k = 1:populationSize
        if(y(k) > maxnew)
            newy(k) = y(k);
            newx(:, :, k) = x(:, :, k);
            val = val + 1;
        else
            newy(k) = nextGenY(k - val);
            newx(:, :, k) = nextGenX(:, :, k - val);
        end
    end
    
    x = newx;
    y = newy;
end

y = f(bin2float(x(:, :, 1)));
end

function res = float2bin(val)
    %Float2bin returns the binary representation of a float number in base 10
    n = 3; %digits before comma
    m = 10; %digits after comma

    res = [val<0 fix(rem(abs(val).*pow2(-(n-1):m),2))];
end

function res = bin2float(val)
    %Bin2float returns the base 10 float representation of a binary number.
    %n = 3;
    %m = 10;
    %pows = pow2(n-1:-1:-m)
    pows = [4.0000    2.0000    1.0000    0.5000    0.2500    0.1250    0.0625    0.0313    0.0156    0.0078    0.0039    0.0020    0.0010];
    s = size(val);
    pows = repmat(pows,s(1),1);
    
    mult = (val(:, 2:s(2)).*pows);
    res = ((-1).^val(:,1)).*sum(mult, 2);
end

function [num1, num2] = crossover(val1, val2)
    %Crossover crosses 2 binary numbers at a random point and returns the result (without changing the inputs)
    digitBefComma = 3;
    digitAftComma = 10;
    
    %Generate random number between 2 and  1 + BinaryNumber.digitBefComma + BinaryNumber.digitAftComma
    uplim = 1 + digitBefComma + digitAftComma;
    digit = ceil(rand(1)*(uplim-2))+1;

    %save the 4 subvectors
    tmp11 = val1(1:digit-1);
    tmp12 = val1(digit:uplim);
    tmp21 = val2(1:digit-1);
    tmp22 = val2(digit:uplim);

    %swap the subvectors
    num1 = [tmp11 tmp22];
    num2 = [tmp21 tmp12];
end

function obj = mutate(obj)
    %GetDigit returns the digit at position n. Returns a new object
    %sign is mutable too

    mutationProb = 0.1;
    digitBefComma = 3;
    digitAftComma = 10;

    %Generate a random number for each digit
    rands = rand(digitBefComma + digitAftComma + 1, 1);

    %mutate every digit where the random number generated is
    %smaller then BinaryNumber.mutationProb
    idxs = find(rands <= mutationProb);
    obj(idxs) = mod(obj(idxs)+1,2);
end

