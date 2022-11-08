function [IMAGE_3D_DATA] = image3Ddata(M)
% function [IMAGE_3D_DATA]=image3Ddata(M)
% ------------------------------------------------------------------------
% 
% This simple function creates a structure array containing coordinate and
% colour data for 3D images. It allows one to use the patch function to
% plot the whole image or a selection of voxels in 3D.
%
% N.B. The function has not been optimised for large images. Large images
% (function has been tested for images under 100x100x100) may produce
% memory problems. 
%
%
% EXAMPLE
%
% M = rand(15,15,15);
% [IMAGE_3D_DATA] = image3D(M);
% 
% Getting faces and vertices for full image:
% voxel_no=1:1:numel(M);
% voxel_face_no=IMAGE_3D_DATA.voxel_patch_face_numbers(voxel_no,:);
% M_faces=IMAGE_3D_DATA.voxel_patch_faces(voxel_face_no,:);
% M_vertices=IMAGE_3D_DATA.corner_coordinates_columns_XYZ;
% 
% Getting faces and vertices for selection of voxels:
% voxel_no2=M>0.95;
% voxel_face_no2=IMAGE_3D_DATA.voxel_patch_face_numbers(voxel_no2,:);
% M_faces2=IMAGE_3D_DATA.voxel_patch_faces(voxel_face_no2,:);
% M_vertices2=IMAGE_3D_DATA.corner_coordinates_columns_XYZ;
% 
% figure;fig=gcf; clf(fig); colordef (fig, 'white'); units=get(fig,'units'); set(fig,'units','normalized','outerposition',[0 0 1 1]); set(fig,'units',units);
% set(fig,'Color',[1 1 1]);
% 
% subplot(1,2,1);
% hp=patch('Faces',M_faces,'Vertices',M_vertices,'EdgeColor','black', 'CData',IMAGE_3D_DATA.voxel_patch_CData(voxel_face_no,:),'FaceColor','flat');
% hold on; view(45,30); axis equal; axis tight; colormap jet; colorbar; caxis([0 1]);
% xlabel('J'); ylabel('I'); zlabel('K');
% title('Full image');
% 
% subplot(1,2,2);
% hp2=patch('Faces',M_faces2,'Vertices',M_vertices2,'EdgeColor','black', 'CData',IMAGE_3D_DATA.voxel_patch_CData(voxel_face_no2,:),'FaceColor','flat');
% hold on; view(45,30); axis equal; axis tight; colormap jet; colorbar; caxis([0 1]); grid on;
% set(hp2,'FaceAlpha',0.8); 
% xlabel('J'); ylabel('I'); zlabel('K');
% title('Selection of voxels');
%
%
% Kevin Mattheus Moerman
% kevinmoerman@hotmail.com
% 11/05/2009
% ------------------------------------------------------------------------


% Copyright (c) 2009, Kevin Mattheus Moerman
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are 
% met:
% 
%     * Redistributions of source code must retain the above copyright 
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright 
%       notice, this list of conditions and the following disclaimer in 
%       the documentation and/or other materials provided with the distribution
%       
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
% POSSIBILITY OF SUCH DAMAGE.

% Setting up meshgrid of voxel centre coordinates
[X,Y,Z] = meshgrid(1:1:size(M,2),1:1:size(M,1),1:1:size(M,3));
IMAGE_3D_DATA.center_coordinates_meshgrid_X=X;
IMAGE_3D_DATA.center_coordinates_meshgrid_Y=Y;
IMAGE_3D_DATA.center_coordinates_meshgrid_Z=Z;

X=X(:); Y=Y(:); Z=Z(:);
IMAGE_3D_DATA.center_coordinates_columns_XYZ=[X Y Z];

% Creating coordinates for voxel corners
[X,Y,Z] = meshgrid(0.5:1:(size(M,2)+0.5),0.5:1:(size(M,1)+0.5),0.5:1:(size(M,3)+0.5));             
             
IMAGE_3D_DATA.corner_coordinates_meshgrid_X=X;
IMAGE_3D_DATA.corner_coordinates_meshgrid_Y=Y;
IMAGE_3D_DATA.corner_coordinates_meshgrid_Z=Z;

X=X(:); Y=Y(:); Z=Z(:);
IMAGE_3D_DATA.corner_coordinates_columns_XYZ=[X Y Z];

clear X Y Z;

% Creating path face and color data
nodes_first_voxel = [    1 ...
                         2 ...
                        ( ( (size(M,1)+1)*(size(M,2)+1) ) +2) ...
                        ( ( (size(M,1)+1)*(size(M,2)+1) ) +1) ...
                        ( 1 + (size(M,1)+1) ) ...
                        ( 2 + (size(M,1)+1) ) ...
                        ( 2 + (size(M,1)+1) + ( (size(M,1)+1)*(size(M,2)+1) ) ) ...
                        ( 1 + (size(M,1)+1) + ( (size(M,1)+1)*(size(M,2)+1) ) )   ];

nodes_first_row_voxels=((0:1:(size(M,1)-1))' * ones(1,8)) + (ones(1,size(M,1))') * nodes_first_voxel;
A = repmat(nodes_first_row_voxels,(size(M,2)),1);
B = repmat(((size(M,1)+1)*(0:(size(M,2)-1))),(size(M,1)),1);
B = reshape(B,1, numel(B))' *ones(1,8);
nodes_first_slice_voxels=A+B;
A = repmat(nodes_first_slice_voxels,(size(M,3)),1);
B = repmat((((size(M,1)+1)*(size(M,2)+1))*(0:(size(M,3)-1))),((size(M,1))*(size(M,2))),1);
B = reshape(B,1,numel(B))'*ones(1,8);

IMAGE_3D_DATA.corner_numbers=A+B;
IMAGE_3D_DATA.voxel_patch_face_numbers=reshape(1:1:(6*numel(M)),6,numel(M))';
IMAGE_3D_DATA.voxel_patch_CData=reshape(((M(:)*ones(1,6)))',(6*numel(M)),1);

face_no=[1 2 3 4;1 2 6 5;2 3 7 6;3 4 8 7;1 4 8 5;5 6 7 8]';
faces=reshape((IMAGE_3D_DATA.corner_numbers(:,face_no))',4,[])';
IMAGE_3D_DATA.voxel_patch_faces=faces;