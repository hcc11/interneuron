function g=WrappedGauss1D(x,mu,sigma,k)
    g=zeros(size(x));
    for j1=-k:k    
        g=g+normpdf(x+j1,mu,sigma);
    end
end