% ****************************************************************************** 
% getVoxelSize - Returns a portion of the spectrum specified by limits (ppm)
%
% [volume, dims] = getVoxelSize(sp)
% returns the volume (in mL) and dimensions (mm) in nominal order

% Created: 11/6/2010 Patrick Bolan
% MBS - Minnesota Breast Spectroscopy package
function [volume, dims] = getVoxelSize(sp)


dims(1) = str2double(sp.voxel.slab(1).thickness);
dims(2) = str2double(sp.voxel.slab(2).thickness);

% Some vendors incorrectly set the other dims 

dims(3) = str2double(sp.voxel.slab(3).thickness);
volume = dims(1)*dims(2)*dims(3) / 1000;