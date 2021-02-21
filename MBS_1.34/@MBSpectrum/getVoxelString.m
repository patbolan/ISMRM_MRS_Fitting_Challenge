% ****************************************************************************** 
% getVoxelString - Returns a string that defines the voxel geometry
%
% str = getVoxelString(sp)
% The DICOM standard voxel coordinates are described awkwardly: a series of
% 3 slabs, that may or may not be intersecting. That's very flexible, but
% very incovenient for comparing 2 voxel geometries when most of the time
% it is simply a cuboid.
% I will use this shorthand showing (position)[dims]{direction cosines}
% (3.2\.3\-23.1)[10.1\12.1\4]{1\0\0\0\1\0}

% Created: 1/26/2012 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function str = getVoxelString(sp)

try 
    slab = sp.voxel.slab;

    % I assume that the position is the same for all
    posstr = slab(1).position;
    
    dimstr = sprintf('[%.1f\\%.1f\\%.1f]', ...
        str2num(slab(1).thickness), ...
        str2num(slab(2).thickness), ...
        str2num(slab(3).thickness) );
     
    orientstr = sprintf('{%s\\%s}', ...
        slab(1).orientation, slab(2).orientation);
    
    str = [dimstr ' ' posstr ' ' orientstr];
    
catch
   str = 'FAILED TO PARSE VOXEL STRING'; 
end

