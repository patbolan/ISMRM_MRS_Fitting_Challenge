function spec3d = extractCSIMatrix(sp)

% Does not support 3D (frames)

nCols = str2num(sp.params.MultivoxelColumns.value{1});
nRows = str2num(sp.params.MultivoxelRows.value{1});

% Make up a 3D struct
spec3d = zeros(nRows, nCols, sp.N);
for mdx=1:sp.M
    [row, col, ~] = sp.indexToCoords(mdx);
    spec3d(row, col, :) = sp.fid(:,mdx);
end



