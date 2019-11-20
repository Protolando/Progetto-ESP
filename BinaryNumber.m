classdef BinaryNumber < handle
    %BinaryNumber Binary representation of a number with up to 13 digits (3 before and 10 after the comma). rep(1) = sign(value)
    %   Class that holds the binary representation of a number and implements functions for binary number manipulation inside a genetic algorithm
    
    properties(Constant)
        digitBefComma = 3;
        digitAftComma = 10;
        mutationProb = 0.1;
    end
    properties
        rep = zeros(1, 1 + BinaryNumber.digitBefComma + BinaryNumber.digitAftComma);
        %+1 for the sign bit
        value = 0;
    end
    
    methods
        function obj = BinaryNumber(number)
            %BinaryNumber Construct an instance of this class from either a float or a bin array
            %   Takes in input a number and stores the binary representation. Returns this.
            if(isscalar(number))
                %If number is scalar it will be treated as a float
                obj.value = number;
                obj.rep = BinaryNumber.Float2bin(number);
            else
                %If number is an array it will be treated as a bit array
                obj.rep = number;
                obj.value = BinaryNumber.Bin2float(number);
            end
        end
        
        function ret = GetDigit(self, n)
            %GetDigit returns the digit at position n
            if (n > BinaryNumber.digitBefComma + BinaryNumber.digitAftComma || n < 1)
                ret = 0;
            else
                ret = self.rep(n);
            end
        end
        
        function obj = Mutate(self)
            %GetDigit returns the digit at position n. Returns a new object
            %sign is mutable too
            
            binNum = self.rep;
            
            %Generate a random number for each digit
            rands = rand(BinaryNumber.digitBefComma + BinaryNumber.digitAftComma + 1, 1);
            
            %mutate every digit where the random number generated is
            %smaller then BinaryNumber.mutationProb
            idxs = find(rands <= BinaryNumber.mutationProb);
            binNum(idxs) = mod(binNum(idxs)+1,2);
            
            obj = BinaryNumber(binNum);
        end
        
    end
    
    methods(Static)
        
        function obj = BinaryNumberArray(size, numbers)
            %BinaryNumber Construct an array of instances of this class from an array of float
            %   takes in input a number of elements and an array, allocates an array of size BinaryNumbers with the values in number
            
            obj(1:size) = BinaryNumber(0);
            for i = 1:size
                obj(i) = BinaryNumber(numbers(i));
            end
        end
        
        function res = Float2bin(val)
            %Float2bin returns the binary representation of a float number in base 10
            n = BinaryNumber.digitBefComma;
            m = BinaryNumber.digitAftComma;
            
            res = [val<0 fix(rem(abs(val).*pow2(-(n-1):m),2))];
        end
        
        function res = Bin2float(val)
            %Bin2float returns the base 10 float representation of a binary number.
            n = BinaryNumber.digitBefComma;
            m = BinaryNumber.digitAftComma;
            
            s = size(val);
            mult = val(2:s(2));
            res = ((-1)^val(1))*mult*pow2(n-1:-1:-m).';
        end
        
        function [num1, num2] = Crossover(val1, val2)
            %Crossover crosses 2 binary numbers at a random point and returns the result (without changing the inputs)
            
            %Generate random number between 2 and  1 + BinaryNumber.digitBefComma + BinaryNumber.digitAftComma
            uplim = 1 + BinaryNumber.digitBefComma + BinaryNumber.digitAftComma;
            digit = ceil(rand(1)*(uplim-2))+1;
            
            %save the 4 subvectors
            tmp11 = val1.rep(1:digit-1);
            tmp12 = val1.rep(digit:uplim);
            tmp21 = val2.rep(1:digit-1);
            tmp22 = val2.rep(digit:uplim);
            
            %swap the subvectors
            num1 = BinaryNumber([tmp11 tmp22]);
            num2 = BinaryNumber([tmp21 tmp12]);
        end
    end
end

