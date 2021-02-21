function [row, col, frame] = indexToCoords(sp, idx)

nCols = str2num(sp.params.MultivoxelColumns.value{1});
nRows = str2num(sp.params.MultivoxelRows.value{1});
numFrames = str2num(sp.params.MultivoxelFrames.value{1});

% First change index for number of frames
if numFrames>1
    sizePerFrame = nRows*nCols;
    frame = floor(idx/sizePerFrame);
    idx = mod(idx, sizePerFrame);
else 
    frame = 1;
end


% This logic does not support 3rd dimension
row = ceil(idx/nCols);
col = mod(idx-1, nCols)+1;


return;
