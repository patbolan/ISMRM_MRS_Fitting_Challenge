% plotSpec Plots the frequency-domain spectrum in the current axes
%
% plotSpec(sp, axis, mode) 
%   Axis = 'h' (Hz, default) or 'p' (ppm)
%   Mode = 0 (real, default), 1 (imag), 2 (abs), other (real+imag).

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotSpec(sp, axis, mode)

if nargin > 1 & axis == 'h'
    freqax = sp.freqHz;
    axlabel = 'Hz';
else
    freqax = sp.freqPPM;
    axlabel = 'ppm';
end

% mode: 0=real, 1=imag, 2=abs, 3=real+imag
if nargin<3
    mode = 0;
end

sp = mean(sp);

switch (mode)
    case 0
        plot(freqax, real(sp.spec), 'b');
    case 1
        plot(freqax, imag(sp.spec), 'b');
    case 2
        plot(freqax, abs(sp.spec), 'b');
    case 3 
        plot(freqax, real(sp.spec), 'b', freqax, imag(sp.spec), 'r');        
    otherwise
        
        % Real + angle
        range = min( max(real(sp.spec)), min(real(sp.spec)) );
        fscale = .5 * range / pi;
        plot(freqax, real(sp.spec), 'b', freqax, angle(sp.spec) .* fscale, '-m');
        legend('real', 'phase');
        
end

set(gca,'Xdir','reverse')
xlabel(axlabel);
zoom on;

