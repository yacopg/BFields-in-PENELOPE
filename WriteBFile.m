% GenBfile.m
%   -script to write a magnetic field map to text file in the format currently expected by
%   the implementation of arbitrary magnetic fields into PENELOPE. 
%   Requires:
%       -BX, BY, BZ: the x, y, and z components of the magnetic field volume. 
%       -x_M, y_M, z_M: vectors defining the x, y, and z spatial coordinates for the
%       magnetic field volume. 

%% Example: A simple magnetic field that will result in a uniformily varying magnetic field from 
% 0 to 10 T over a distance of 50 cm.
%coordinate grid.
z_M = [0,0.5];
y_M = [-0.15,0,0.15];
x_M = [-0.15,0,0.15];

sz = [length(y_M), length(x_M), length(z_M)];

BX = zeros(sz); %x-component of B field
BY = zeros(sz);
BZ = zeros(sz);

BZ(:,:,1) = 0;
BZ(:,:,2) = 10;

%% Write out the magnetic field map. 
BOutFile = "test.txt";

fID = fopen(BOutFile,'w');
fprintf(fID,'%6.8f %6.8f %6.8f\n',x_M(1),y_M(1),z_M(1)); % first coordinates
fprintf(fID,'%6.8f %6.8f %6.8f\n',x_M(2) - x_M(1),y_M(2) - y_M(1), ...
                                  z_M(2) - z_M(1)); % step lengths
fprintf(fID,'%6i %6i %6i\n',length(x_M),length(y_M),length(z_M)); % N
fprintf(fID,'%6.8f\n',BX);
fprintf(fID,'%6.8f\n',BY);
fprintf(fID,'%6.8f\n',BZ);
fclose(fID);

