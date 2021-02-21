% ******************************************************************************
%             MBS - Minnesota Breast Spectroscopy analysis package
%               Developed by Patrick Bolan and Michael Garwood
% ******************************************************************************
% FUNCTION: mbsSpectrum - plotSpecArrayMovie
% AUTHOR: pjb
% CREATED: 8/7/2002
% DESCRIPTION: Plots the summed spec
% ARGUMENTS: 
%   sp - mbsSpectrum
%   axis - ('h' for Hz, 'ppm'), defaults to ppm
%   mode - 0=real, 1=imag, 2=abs, 3=real+imag (defaults to 3)
%   xlim - defaults to Matlab's plot default for the first scan
%   ylim - defaults to Matlab's plot default for the first scan
% RETURNS: none
% MODIFICATIONS:
% ******************************************************************************

% plotSpecArrayMovie - Shows a movie of all the spectra in the array
%
% plotSpecArrayMovie(sp, axis, mode, xlim, ylim, delay) 
%   axis - ('h' for Hz, 'ppm'), defaults to ppm
%   mode - 0=real, 1=imag, 2=abs, 3=real+imag (defaults to 3)
%   xlim - defaults to Matlab's plot default for the first scan
%   ylim - defaults to Matlab's plot default for the first scan

% Created: 8/7/2002 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function plotSpecArrayMovie(sp, axis, mode, xlim, ylim, delay)

if nargin > 1 & axis == 'h'
    freqax = sp.freqHz;
    axlabel = 'Hz';
else
    freqax = sp.freqPPM;
    axlabel = 'ppm';
end

% mode: 0=real, 1=imag, 2=abs, 3=real+imag
if nargin<3
    mode = 3;
end

    
spone = extract(sp, 1);

switch (mode)
    case 0
        plot(freqax, real(spone.spec), 'b');
    case 1
        plot(freqax, imag(spone.spec), 'b');
    case 2
        plot(freqax, abs(spone.spec), 'b');

    otherwise
        plot(freqax, real(spone.spec), 'b', freqax, imag(spone.spec), 'g');
end

set(gca,'Xdir','reverse')
xlabel(axlabel);

if (nargin<4)
    bAutoscale = 1;
else
    bAutoscale = 0;
    set(gca, 'Xlim', xlim);
    set(gca, 'Ylim', ylim);
end

time_delay = 0.5;
if (nargin > 5 )
    % Fix the delay to 0.5s
    time_delay = delay;
end
    

for idx = 2:sp.M
    %dummy = input(sprintf('Showing %d/%d (press return)', idx-1, numspec));
    strMsg = sprintf('Spectrum %d/%d', idx-1, sp.M);
    title(strMsg)
    pause(time_delay);

    spone = extract(sp, idx);

    switch (mode)
        case 0
            plot(freqax, real(spone.spec), 'b');
        case 1
            plot(freqax, imag(spone.spec), 'b');
        case 2
            plot(freqax, abs(spone.spec), 'b');

        otherwise
            plot(freqax, real(spone.spec), 'b', freqax, imag(spone.spec), 'g');
    end

    set(gca,'Xdir','reverse')
    if (bAutoscale == 0)
        set(gca, 'Xlim', xlim);
        set(gca, 'Ylim', ylim);
    end
    xlabel(axlabel);

end
