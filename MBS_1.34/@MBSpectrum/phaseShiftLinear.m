% phaseSpecLinear - perform a linear (first order) phase shift
%
% sp = phaseSpecLinear(sp, angDegrees, lpDegrees) - perform
%   angRadians - this operates in one of two modes. If angRadians is scalar, 
%   it phases all spectra the same. If its an array, it better be equal to the 
%   number of spectra - then each gets different phases. In radians
%   lpRadians is interpreted as radians per ppm.

% Created: 2/6/2006 Patrick Bolan,
% MBS - Minnesota Breast Spectroscopy package
function sp = phaseShiftLinear(sp, angRadians, lpRadians)

dim = max(size(angRadians));
if dim ==1
    angRadians = (1:sp.M).* 0 + angRadians;
    lpRadians = (1:sp.M).* 0 + lpRadians;
elseif max(size(angRadians)) ~= sp.N
    error('Dimension of angRadians is not equal to the array size M');  
end


% This introduces some artifacts at the start and end of the FID. Need to
% revisit the linear phase algorithm
for idx = 1:sp.M

	% Performance - zero order is faster
	if (lpRadians(idx) == 0) 
		sp.spec(:,idx) = sp.spec(:,idx) .* exp(-i*angRadians(idx));
	else
		% Calculate the phase function
		clear phs;
        
        phs = angRadians(idx) + (((sp.N/2)-(1:sp.N))/sp.N) * (lpRadians(idx)*sp.swppm);
        
        % This old way was very slow
% 		for jdx = 1:sp.N
% 			phs(jdx) = angRadians(idx) + ...
%                 (((sp.N/2)-jdx)/sp.N) * (lpRadians(idx)*sp.swppm);
% 		end
		phs = phs';

		% Apply it
		sp.spec(:,idx) = sp.spec(:,idx) .* exp(-i.*phs);
	end
end














