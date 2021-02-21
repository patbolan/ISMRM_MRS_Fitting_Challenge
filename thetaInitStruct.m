function theta = thetaInitStruct(Nregions, Nelements)

% globals
theta.global(1).name = 'phs'; % degrees
theta.global(1).value = 0.0;
theta.global(2).name = 'llw'; % Hz
theta.global(2).value = 0.0;
theta.global(3).name = 'glw'; % Hz
theta.global(3).value = 0.0;
theta.global(4).name = 'delf'; % Hz
theta.global(4).value = 0.0;


for rdx=1:Nregions
    for edx=1:Nelements
        theta.region(rdx).element(edx).label = '';
        theta.region(rdx).element(edx).amp = 0.0;
        theta.region(rdx).element(edx).phs = 0.0; % degrees
        theta.region(rdx).element(edx).llw = 0.0; % Hz
        theta.region(rdx).element(edx).glw = 0.0; % Hz
        theta.region(rdx).element(edx).delf = 0.0; % Hz        
    end
end

