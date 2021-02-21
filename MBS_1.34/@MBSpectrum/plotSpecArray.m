% plotSpecArray Plots an array of frequency-domain spectra in the current axes
%
% plotSpecArray(sp, axis, hoffset, voffset, mode)) 
%   Axis = 'h' (Hz, default) or 'p' (ppm)
%   Mode = 0 (real, default), 1 (imag), 2 (abs), other (real+imag).

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotSpecArray(sp, axis, hoffset, voffset, mode)
cla
if nargin > 1 & axis == 'h'
    freqax = sp.freqHz;
else
    freqax = sp.freqPPM;
end

if nargin <= 2
    hoffset = 0;
    voffset = 1;
end

% mode: 0=real, 1=imag, 2=abs
if nargin<5
    mode = 0;
end

spave = mean(sp);
offset = max(abs(spave.spec)) * voffset/100;
% Changing to absolute offset
offset= voffset;

colorOrder = get(gca, 'ColorOrder');

% Hey, I think I could optimize this but making 1 plot call, but
% prepaing a matrix of values
hold on
for idx = 1:sp.M
    switch (mode)
        case 0
            signal = real(sp.spec(:,idx));

        case 1
            signal = imag(sp.spec(:,idx));

        otherwise
            signal = abs(sp.spec(:,idx));
    end
    
    coff =  mod(idx-1,size(colorOrder,1))+1;
    plot(freqax - (idx-1)*hoffset, ...
        signal+(idx-1)*offset, 'Color', colorOrder(coff, :));
end
hold off
set(gca,'Xdir','reverse')
