% plotFidArray Plots array of the time-domain spectra in the current axis
%
% plotFidArray(sp, hoffset, voffset, mode)
%   mode: 0=real, 1=imag, 2=abs (default), 3=real+imag

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotFidArray(sp, hoffset, voffset, mode)
clf
timeax = sp.time;

if nargin < 2
    hoffset = 0;
    voffset = .01;
end

% mode: 0=real, 1=imag, 2=abs
if nargin<4
    mode = 0;
end

fidave = mean(sp);
offset = max(abs(fidave.fid)) * voffset;

colorOrder = get(gca, 'ColorOrder');

hold on
for idx = 1:sp.M
    switch (mode)
        case 0
            signal = real(sp.fid(:,idx));
            
        case 1
            signal = imag(sp.fid(:,idx));
            
        otherwise
            signal = abs(sp.fid(:,idx));
    end
    
    coff =  mod(idx-1,size(colorOrder,1))+1;
    plot(timeax + (idx-1)*hoffset, ...
        signal+(idx-1)*offset,'Color', colorOrder(coff, :));
end
hold off
if offset>0
    set(gca,'Ylim', [-5*offset 150*offset])
end
