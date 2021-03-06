function h = draw(G,varargin)    

varlist = varargin;

patch_option = find(cellfun(@(x) isstruct(x), varlist));
if ~isempty(patch_option)
    options = varlist{patch_option};
    varlist(patch_option) = [];
end


if ~isempty(find(cellfun(@(x) strcmpi(x,'flat'), varlist), 1))
    MeshType = 'flat';
    varlist(cellfun(@(x) strcmpi(x,'flat'), varlist)) = [];
else
    MeshType = 'original';
end

if ~isempty(varlist)
    FeatureType = varlist;
end

F = G.F;
V = G.V;
if strcmpi(MeshType,'flat')
    V = G.Aux.UniformizationV;
end

if isempty(F) || isempty(V)
    disp('Empty vertices or faces!');
    return;
end
if ~iscell(F) %convert to cell array
    F=mat2cell(F',ones(1,size(F,2)),3);
end

face_size = cellfun(@numel, F); % get face size per face
minmax_face_size = [min( face_size ),max(face_size )];

if minmax_face_size(1)<3
    warning('Non polygonal faces: %s', num2str(find(fdeg<3)'));
    minmax_face_size(1) = 3;
end

h = zeros( minmax_face_size (2), 1 );

if size(F{1},1)~=1
    F = cellfun( @(x) x', F, 'UniformOutput', 0);
end

for i=minmax_face_size(1):minmax_face_size(2)
    h(i) = patch('Vertices', V', 'Faces', cell2mat(F(face_size==i)),'FaceColor',[0.5,0.5,0.5],'CDataMapping','direct');
end

h = h( h>0 );
if exist('options', 'var')==1 && ~isempty(options)
    set(h, options);
end

dim = size(V, 1);
if dim==3 && max(abs(V(3,:))) < 1e-8
    dim = 2;
end

if dim<3
    view(2);
else
    view(3);
end

grid off;
axis off;
axis tight;
axis equal; 

if exist('FeatureType','var')
    hold on;
    for j=1:length(FeatureType)
        switch lower(FeatureType{j})
            case 'confmax'
                scatter3(V(1,G.Aux.ConfMaxInds),V(2,G.Aux.ConfMaxInds),V(3,G.Aux.ConfMaxInds),30,'g','filled');
            case 'gaussmax'
                scatter3(V(1,G.Aux.GaussMaxInds),V(2,G.Aux.GaussMaxInds),V(3,G.Aux.GaussMaxInds),30,'r','filled');
            case 'gaussmin'
                scatter3(V(1,G.Aux.GaussMinInds),V(2,G.Aux.GaussMinInds),V(3,G.Aux.GaussMinInds),30,'b','filled');
            case 'admax'
                scatter3(V(1,G.Aux.ADMaxInds),V(2,G.Aux.ADMaxInds),V(3,G.Aux.ADMaxInds),30,'y','filled');
            otherwise
                disp('Un-recognizable Feature Type');
        end
    end
end

cameratoolbar;
cameratoolbar('SetCoordSys', 'none');

