% Takes the theta struct and pretty-prints it out.
function thetaPrint(theta)

nGlobals = size(theta.global,2);
nRegions = size(theta.region,2);
nElements = size(theta.region(1).element, 2);

fprintf('Theta structure with %d globals, %d regions, %d elements per region \n',...
    nGlobals, nRegions, nElements);

% Generic globals in native units
for gdx=1:nGlobals
    fprintf(' \tGlobal %s %f\n', theta.global(gdx).name, theta.global(gdx).value);
end


for rdx=1:nRegions
    fprintf('Region %d\n', rdx)
    
    for edx = 1:nElements
        
        fprintf('\t%d (%s): amp %.3f, \tphs %.3f deg, \tllw %.3f Hz, \tglw %.3f Hz, \tdelf %.3f Hz\n', ...
            edx, ...
            theta.region(rdx).element(edx).label, ...
            theta.region(rdx).element(edx).amp, ...
            theta.region(rdx).element(edx).phs, ...
            theta.region(rdx).element(edx).llw, ...
            theta.region(rdx).element(edx).glw, ...
            theta.region(rdx).element(edx).delf)
        
    end
end





