function fij=fitness(zi,z0,fmax,k)
% compute the fitness of j-th species at the elevation z = zi
% Logistic growth function
fij=(fmax)./(1+exp(-k.*(zi-z0)));

end