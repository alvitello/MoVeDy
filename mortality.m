function muij=mortality(zi,mu0,z0,fmax,k)
% compute the mortality rate of j-th species at the elevation z = zi

muij=mu0./(fitness(zi,z0,fmax,k)./fmax);
end