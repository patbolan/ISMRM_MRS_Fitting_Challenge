function thetaArr = thetaStructToArray(theta)

thetaArr = [];
cnt = 0; % Current index. Treating like queue

% first are the glolbals
nGlobals = size(theta.global,2);
for idx=1:nGlobals
    thetaArr(cnt+1) = theta.global(idx).value;
    cnt = cnt+1;
end


nRegions = size(theta.region,2);
nElements = size(theta.region(1).element, 2);

for rdx = 1:nRegions
    for edx = 1:nElements
        thetaArr(cnt+1) = theta.region(rdx).element(edx).amp;
        thetaArr(cnt+2) = theta.region(rdx).element(edx).phs;
        thetaArr(cnt+3) = theta.region(rdx).element(edx).llw;
        thetaArr(cnt+4) = theta.region(rdx).element(edx).glw;
        thetaArr(cnt+5) = theta.region(rdx).element(edx).delf;
        cnt = cnt + 5;
    end
end



