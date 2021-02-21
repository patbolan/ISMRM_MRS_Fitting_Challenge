classdef MBSpectrum
    % A class that represents an array of 1D spectra.
    % This is designed to be compatible with standard DICOM MRS
    % Please see the documentation in MBS/doc/home.html

    properties
        % The fid is the fundamental data. It is an NxM array of the
        % complex, double, time-domain data, where N is the number of
        % points and M is the size of the array of spectra (ie number of
        % traces)
        fid = [];

        % These are the few pieces of essential header information that
        %   are stored that define the physical properties of the spectrum.

        frq = 100;      % transmitter/receiver FReQuency in MHz
        ppmref = 0;     % Frequency of center, in ppm
        sw = 1000;      % Full SpectralWidth of spectrum in Hz


        % The header information that is read in is stored in 2 structures,
        % the params and the header. The header is a very loose and
        % undefined set of tag-value pairs. This is so that any file reader
        % can put their data into this structure. These values are not
        % written out to any of the standard file writers.
        header = [];    % This is a very loose set of tag-value pairs

        % The params data structure has a specific format so that it can
        % written out in a structured way, either for DICOM MRS or the XML
        % format. 
        % A parameter has the following fields:
        %  name - this is the field name, in matlab
        %  value - always stored as a string, but may be interpreted as a
        %   number. In memory, these are Matlab Cells, with size 1x1 for a
        %   singular value, and 1xM for an arrayed value
        %  unit - an optional, singular string
        %  type - an optional, singluar string with describes the name.
        %    An example is name="UniqueIdentifier", which can have type =
        %    "GE RUN NUMBER" or "DICOMUID"
        %  arrayed - an optional, singular string, true or false. If
        %   arrayed==true then the size of the value array must be equal to M
        % 
        % params is an awkward construct but was necessary to handle writing 
        % of arrayed parameters.
        % Note that arrayed parameters are preserved when two objects are
        % appended or extracted IFF the "a" parameter is arrayed. See
        % append for details
        
        params = [];    

        % The voxel localization is a special structure
        voxel = [];
    end

    % These are read-only properites that are calculated on demand. In
    % previous versions of MBS these values were all saved and they were
    % kept in sync. This approach is a bit slower but cleaner.
    properties (Dependent = true)
        spec    % A NxM array of the complex, double, freq-domain data
        time    % An array of length N holding the time values in s
        freqPPM    % An array of length N holding the freq values in ppm
        freqHz  % An array of length N holding the freq values in Hz

        dt      % Temporal sampling period, s
        at      % Total time of aquisition, s
        dfrq    % Spectral resolution, Hz
        N   	% Number of points in fid/spec
        M   	% Size of array (aka number of traces)
        swppm   % Spectral width in ppm
    end


    methods
        % Constructor
        function obj = MBSpectrum(s)
            if nargin<1
                % Default constructor
            elseif isa(s, 'MBSpectrum');
                % TODO: copy all internal values
            else
                error('cannot create an object with a %s', class(s));
            end
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Getting and setting the spectral data requires FFT/iFFT from the
        % stored FID
        function spec = get.spec(sp)
            % The spectral data is calculated from the fid on demand
            spec = sp.fid .* 0;

            % This is the fast, matrix form of fft + fftshift
            %spec = fft(sp.fid,[],1);
            %shiftidx = cat(2,[sp.N/2+1:sp.N],[1:sp.N/2]);
            %spec(1:sp.N,:) = spec(shiftidx,:);

            % This is the standard implementation
            %spec = fftshift(fft(sp.fid), 1);
            
            % This scales the first point by 1/2 to correct the DFT. See
            % Zhu et al, JMR A 105, 1993
            tmp = sp.fid;
            tmp(1,:) = tmp(1, :) * 0.5;
            spec = fftshift(fft(tmp), 1);            
            
            % Normalize the frequency spectrum by the number of points
            %   so that the units of the time and frequency domain are
            %   equal.
            spec = spec ./ sqrt(sp.N);
        end

        function sp = set.spec(sp, spec)
            sp.fid = spec .* 0;
            
            % In prior versions I had some trouble with fftshift, but I
            % don't see this problem in Matlab 7.6 or greater
            sp.fid = ifft(fftshift(spec, 1));
            
            % And correct for the 0.5 scaling, as in get.spec
            sp.fid(1, :) = sp.fid(1, :) .* 2;

            % Adjust the magnitude here to keep the fid and spec in
            %   same scale.
            sp.fid = sp.fid .* sqrt(sp.N);
            
         end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dependent (calculated) fields that are read-only
        function time = get.time(sp)
            % Calculate the time axes, s
            %timeold = (0:sp.at/(sp.N-1):sp.at)';
            time = (0:sp.N-1)' .* sp.dt;
        end

        function freq = get.freqHz(sp)
            % Calculate the frequency axis, Hz

            freq = (sp.sw/2:-sp.dfrq:-sp.sw/2)';
            %freq = freq ./ sp.frq;
            freq = freq + (sp.ppmref*sp.frq);
            %freq = freq + sp.ppmref;
        end
        
        function freq = get.freqPPM(sp)
            % Calculate the frequency axis, ppm

            freq = (sp.sw/2:-sp.dfrq:-sp.sw/2)';
            freq = freq ./ sp.frq;
            freq = freq + sp.ppmref;
        end        

        function pts = get.N(sp)
            % Returns the number of points in each spectrum
            pts = size(sp.fid,1);
        end

        function pts = get.M(sp)
            % Returns the number of spectra in the array
            pts = size(sp.fid,2);
        end

        function at = get.at(sp)
            at = sp.dt * sp.N;
        end

        function dt = get.dt(sp)
            % Returns the temporal resolution, in seconds
            dt = 1/sp.sw;
        end

        function dfrq = get.dfrq(sp)
            % Returns the frequency resolution, in Hz
            if sp.N>0
                dfrq = sp.sw/(sp.N-1);
            else
                dfrq = 0;
            end
        end

        function swppm = get.swppm(sp)
           % Returns the spectral width in ppm units
           swppm = sp.sw / sp.frq;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Array access
        function sp = extract(sp, indices)
            sp.fid = sp.fid(:,indices);
                      
            % Handle arrayed parameters
            oldparams = sp.params;
            arrayedParameters = sp.getArrayedParameters();
            for jdx=1:size(arrayedParameters,2)
                name = arrayedParameters{jdx};
                sp.params.(name).value = { oldparams.(name).value{indices} };
            end
        end

        function a = append(a, b);
            % Appends spectrum 'b' to the end of sp, so that the new array
            %    size is a.M + b.M. 
            % This only works if the array size is correct.
            % The header and parameters from 'a' are dominant. 
            % If any of the 'a' parameters are arrayed, then the appending
            %  will append the b parameter to that array.
            % Note that if other header info is different (frq, etc) the
            %   results will not make any sense. User beware.
            if isa(b, 'MBSpectrum') && (a.N == b.N)
                a.fid(:, a.M + (1:b.M)) = b.fid;
                
                % Find any 'a' parameters that are arrayed
                names = fields(a.params);
                for idx = 1:size(names,1)
                    
                    % Some fields are not arrayable
                    if (    (strcmpi(names{idx}, 'SeriesDescription') == 1) || ...
                            (strcmpi(names{idx}, 'SeriesDate') == 1) || ...
                            (strcmpi(names{idx}, 'SeriesTime') == 1) || ...
                            (strcmpi(names{idx}, 'SourceFile') == 1) )
                        % Skip these - should not be arrayed
                        1;
                        
                    elseif (isfield(a.params.(names{idx}), 'arrayed')) && ...
                            (strcmpi(a.params.(names{idx}).arrayed, 'true') == 1)
                        
                        % Add the b values to that array
                        if ( isfield(b.params, names{idx}) )
                            a.params.(names{idx}).value = ...
                                {a.params.(names{idx}).value{:} b.params.(names{idx}).value{:}};
                        else
                            % IT wasn't in B! Just add a null string
                            a.params.(names{idx}).value = ...
                                {a.params.(names{idx}).value{:} ''};
                            
                        end
                        
                    elseif ( isfield(b.params, names{idx})  && ...
                            strcmp(a.params.(names{idx}).value{1}, b.params.(names{idx}).value{1}) ~=1 )
                        % Alternately, if the b value is different and not
                        % null, make this arrayed                      
                        
                        % Convert parameter in a to array.
                        
                        % First, if the size of values is not really
                        % arrayed, then array it
                        vals = a.params.(names{idx}).value{:};
                        % TODO!
                        
                        % Now amke it arrayed
                        a.params.(names{idx}).arrayed = 'true';
                        a.params.(names{idx}).value = ...
                                {a.params.(names{idx}).value{:} b.params.(names{idx}).value{:}};
                        
                    end
                end
            end
        end
        
        function arrayednames = getArrayedParameters(a)
            if isempty(a.params)
                arrayednames = [];
            else
                names = fields(a.params);
                arrayednames = {};
                for idx = 1:size(names,1)
                    
                    if (isfield(a.params.(names{idx}), 'arrayed')) && ...
                            (strcmpi(a.params.(names{idx}).arrayed, 'true') == 1)
                        arrayednames = {arrayednames{:} names{idx}};
                    end
                end
            end
        end
        
        function sp = setParameter(sp, name, value, arrayed, unit, type)
            % Helper to make it easier to correctly set parameters.
            if nargin<6; type = []; end;
            if nargin<5; unit = []; end;
            if nargin<4; arrayed = []; end;
            
            if nargout<1
               error('You must assign an output value for setParameter to work'); 
            end
            
            % Doesn't matter if the parameter already exists
            
            if isnumeric(value)
                if (max(size(value)) > 1)
                    % Arrayed number
                    arrayed = 'true';
                    cellval = cell(1,max(size(value)));
                    for idx=1:max(size(value))  
                        cellval{1,idx} = num2str(value(idx));
                    end
                else
                    % singular number
                    cellval = {num2str(value)};
                end
            elseif ischar(value)
                cellval = {value};
            elseif iscell(value)
                cellval = value;
            end
            
            if ~iscell(cellval); error('Cannot convert to cell.'); end
            
            % set it
            sp.params.(name).value = cellval;
            if ~isempty(arrayed); sp.params.(name).arrayed = arrayed; end;
            if ~isempty(unit); sp.params.(name).unit = unit; end;
            if ~isempty(type); sp.params.(name).type = type; end;
            
        end
        
        % Alternate syntax. Actually longer.
        function [value, arrayed, unit, type] = getParameter(sp, name)
            value = sp.params.(name).value;
            if isfield(sp.params.(name), 'arrayed');
                arrayed = sp.params.(name).arrayed;
            else
                arrayed = [];
            end
            if isfield(sp.params.(name), 'unit');
                unit = sp.params.(name).unit;
            else
                unit = [];
            end
            if isfield(sp.params.(name), 'type');
                type = sp.params.(name).type;
            else
                type = [];
            end
        end
                    
        % Not implementing subsref. It would lead to a nicer subscripting
        % syntax, but then you need to have a switch/case for the different
        % methods and properties (see DocPolynum for example), which would
        % be a pain to maintain.
        %function b = subsref(a, s)


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Several trivial processing methods. More complicated methods are
        % kept in other files

        function sp = mean(sp)
            % Averages the spectra in the array
            sp.fid = mean(sp.fid, 2);
        end

        function sp = sum(sp)
            % Sums the spectra in the array
            sp.fid = sum(sp.fid, 2);
        end
        
        function sp = times(sp, b)
            % Scales the full spectrum/fid by a scalar value. Works in either
            % if b=scalar, value is multiplied to all points in the fid/spec
            % if b=spectrum with 1 spectrum, it is multiplied
            % point-by-point
            % if b=spectrum with mulitple spectra, and it has the same
            % array size as a, they are mulitplied spectrum-by-spectrum and
            % point-by-point

            if (isa(b, 'MBSpectrum'))
                if b.M == 1
                    % Multiply b's spectrum to each of sp's
                    fid = sp.fid;
                    for idx=1:sp.M
                        fid(:,idx) = fid(:,idx) .* b.fid;
                    end
                    sp.fid = fid;

                elseif (sp.M == b.M)
                    % Add each array element pairwise
                    fid = sp.fid;
                    for idx=1:sp.M
                        fid(:,idx) = fid(:,idx) .* b.fid(:,idx);
                    end
                    sp.fid = fid;

                else
                    error('Cannot multiply spectra with unequal array sizes');
                end
            elseif isscalar(b)
                sp.fid = sp.fid .* b;
            end

        end


        % mtimes is not implemented
        function sp = mtimes(sp, b)
            error('mtimes is not implemented. Use ".*" instead');
        end

        % Scalar division is reasonable.
        function sp = rdivide(sp, b);
            if (isa(b,'MBSpectrum'));
                error('Cannot divide by a spectrum, only scalars');
            elseif isscalar(b) && (b~=0)
                sp.fid = sp.fid ./ b;
            else
                error('Unsupported division');
            end

        end

        % This is used to scale the first point to correct for baseline
        function sp = scaleFirstTimepoint(sp, b)
            if (isa(b,'MBSpectrum'));
                error('Cannot scake by a spectrum, only scalars');
            else
                sp.fid(1,:) = sp.fid(1,:) .* b;
            end         
        end
        
        
        function sp = plus(sp, b)
            % Addition. The first argument must be a spectrum.
            % if b=scalar, value is added to all points in the fid/spec
            % if b=spectrum with 1 spectrum, it is added to all
            % if b=spectrum with mulitple spectra, and it has the same
            % array size as a, they are added spectrum-by-spectrum
            if (isa(b, 'MBSpectrum'))
                if b.M == 1
                    % Add b's spectrum to each of sp's
                    fid = sp.fid;
                    for idx=1:sp.M
                        fid(:,idx) = fid(:,idx) + b.fid;
                    end
                    sp.fid = fid;

                elseif (sp.M == b.M)
                    % Add each array element pairwise
                    fid = sp.fid;
                    for idx=1:sp.M
                        fid(:,idx) = fid(:,idx) + b.fid(:,idx);
                    end
                    sp.fid = fid;

                else
                    error('Cannot add spectra with unequal array sizes');
                end
            elseif isscalar(b)
                sp.fid = sp.fid + b;
            end

        end

        % Subtraction is just inverse of addition
        function sp = minus(sp, b)
            sp = plus(sp, b.*-1);
        end

        function sp = phaseShift(sp, phs)
            % Phases the spectra by phs, specified in radians. If phs is a
            % scalar, all spectra will be phased identically. If phs is an
            % array of the same size as the spectrum array, each spectrum
            % will be phased by the distinct phase values

            dim = max(size(phs));
            if dim == 1
                phs = (1:sp.M) .* 0 + phs;
            elseif max(size(phs)) ~= sp.M
                error('Dimension of phases is not equal to array size');
            end

            % This version is very slow! Need to make a matrix version
            phasor = sp.spec .* 0;
            for idx = 1:sp.M
               phasor(:, idx) = exp(i*phs(idx));
            end
            sp.spec = sp.spec .* phasor;

        end
    end

end