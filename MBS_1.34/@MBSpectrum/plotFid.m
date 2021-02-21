% plotFid Plots the time-domain spectrum in the current axis
%
% plotFid(sp, mode) 
%   mode: 0=real, 1=imag, 2=abs (default), 3=real+imag

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotFid(sp, mode)

if nargin<2
    mode = 2;
end

sp = mean(sp);

switch (mode)
    case 0
        plot(sp.time, real(sp.fid), 'b');
        ylabel('real');
    case 1
        plot(sp.time, imag (sp.fid), 'b');
        ylabel('imaginary');

    case 2
        plot(sp.time, abs(sp.fid), 'b');
        ylabel('magnitude');

    case 3
        plot(sp.time, real(sp.fid), 'b', sp.time, imag(sp.fid), 'g');
        legend('real', 'imag');

    case 4
        plot(sp.time, phase(sp.fid), '-m');
        ylabel('phase (rad)');
end


xlabel('Time (ms)');
zoom on;

