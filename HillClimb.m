function x = HillClimb(f, x0, ndim)
   %Initialization
   stepSize = 0.1;
   acceleration = 1.2;
   niters = 10*ndim; %10 iteration x number of dimentions
   xn = x0;
   stopcond = 0.0001; %stop if diff between 2 iteration is less then this
   
   diff = 1;
   i=1;
   
   while (i <= niters && diff > stopcond)
       %Move in a random direction for each dimension
       xn1 = xn + stepSize*(-1).^round(rand(ndim, 1));
       
       %if the new position is better then the old then keep it
       if(f(xn1) > f(xn))
           diff = mean(abs(xn1-xn));
           xn = xn1;
       end
       
       %decrease step size
       stepSize = stepSize/acceleration;
       i = i + 1;
   end
   
   x = xn;
end