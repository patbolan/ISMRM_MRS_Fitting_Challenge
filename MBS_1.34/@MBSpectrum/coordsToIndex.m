function idx = coordsToIndex(sp, row, col, frame)
nCols = str2num(sp.params.MultivoxelColumns.value{1});
nRows = str2num(sp.params.MultivoxelRows.value{1});
%numFrames = str2num(sp.params.MultivoxelFrames.value{1});
idx = (row-1)*nCols + col;

% Adjust for 3D
if (nargin>3)
    frameoffset = (frame-1) * (nRows*nCols);
    idx = idx +frameoffset;
end
return;
