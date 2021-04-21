function [outputArg1] = sumbessel3(x,k,t)
%%%determines the bessel sums used in the analytic solutions
lambdak = besselzero(1,k);%finds the zero of the bessel first order function to be inputted into the zero order bessel function below
tobesummed = zeros(1,length(lambdak));

for i = 1:k  
    tobesummed(i) = besselj(0,x*lambdak(i))*exp(-(t*lambdak(i)^2))/lambdak(i)^2/besselj(0,lambdak(i));  
end

outputArg1 = sum(tobesummed);

end