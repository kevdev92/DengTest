function x = Cintfun(r,t,k)
%%% determines the integral of the concentration times the radius which we
%%% have determined to be given by x below
lambdak = besselzero(1,k); %finds the zero of the bessel first order function to be inputted into the zero order bessel function below
tobesummed = zeros(1,length(lambdak));

for i = 1:k   
    tobesummed(i) = besselj(1,r*lambdak(i))*exp(-(t*lambdak(i)^2))/lambdak(i)^3/besselj(0,lambdak(i));    
end

x = 0.125*(r)^4 + (t-0.125)*(r)^2 - 2*r*sum(tobesummed);

end

