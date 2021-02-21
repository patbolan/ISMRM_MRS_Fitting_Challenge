function theta = thetaArrayToStruct(thetaArr, ref)

theta = ref;
cnt = 0;

% first are the glolbals
nGlobals = size(theta.global,2);
for idx=1:nGlobals
    theta.global(idx).value = thetaArr(cnt+1);
    cnt = cnt+1;
end


nRegions = size(theta.region,2);
nElements = size(theta.region(1).element, 2);

for rdx = 1:nRegions
    for edx = 1:nElements
        theta.region(rdx).element(edx).amp = thetaArr(cnt+1);
        theta.region(rdx).element(edx).phs = thetaArr(cnt+2);
        theta.region(rdx).element(edx).llw = thetaArr(cnt+3);
        theta.region(rdx).element(edx).glw = thetaArr(cnt+4);
        theta.region(rdx).element(edx).delf = thetaArr(cnt+5);
        
        cnt = cnt + 5;
    end
end



