function [vertex,cells,idx,DM,LM] = read_vtk(filename, verbose)

% read_vtk - read data from VTK file.
%
%   [vertex,face] = read_vtk(filename, verbose);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
%   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
%
%   Copyright (c) Mario Richtsfeld

if nargin<2
    verbose = 1;
end

fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(3:5), 'vtk')
    error('The file is not a valid VTK one.');    
end

%%% read header %%%
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);
nvert = sscanf(str,'%*s %d %*s', 1);

% read vertices
[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;

str = fgets(fid);


% read cells
str = fgets(fid);
ncells = sscanf(str,'%*s %d %*s', 1);

[B,cnt] = fscanf(fid,'%f %f %f %f %f', 5*ncells);
if cnt~=5*ncells
    warning('Problem in reading vertices.');
end
B = reshape(B, 5, cnt/5);
cells = B;

str = fgets(fid);
str = fgets(fid);


%Skip
[~,~] = fscanf(fid,'%f %f %f %f %f', 5*ncells);
str = fgets(fid);
str = fgets(fid);
str = fgets(fid);

%IDs
[C,cnt] = fscanf(fid,'%f', nvert);
idx = C;


%Get DarkMem
str = fgets(fid);
str = fgets(fid);
[C,cnt] = fscanf(fid,'%f', nvert);
DM = C;

%Get LitMem
str = fgets(fid);
str = fgets(fid);
[C,cnt] = fscanf(fid,'%f', nvert);
LM = C;


fclose(fid);

return