function [] = ViewObLmk(Names,options)

sample_path = options.sample_path;
DisplayLayout = options.DisplayLayout;
GroupSize = length(Names);
mesh_list = cell(size(Names));
R = options.R;
linkCamera = getoptions(options,'linkCamera','on');
DisplayOrient = getoptions(options,'DisplayOrient','Horizontal');

switch DisplayOrient
    case 'Vertical'
        DisplayOrder = reshape(1:DisplayLayout(1)*DisplayLayout(2), DisplayLayout(2), DisplayLayout(1));
        DisplayOrder = DisplayOrder';
        DisplayOrder = DisplayOrder(:);
    case 'Horizontal'
        DisplayOrder = 1:DisplayLayout(1)*DisplayLayout(2);
end

for i=1:GroupSize
    GM = load([sample_path Names{i} '.mat']);
    GM = GM.G;
    % align every tooth to the first one on the list
    if (i==1)
        mesh_list{i} = GM;
    else
        GM.V = R{1,i}*GM.V;
        mesh_list{i} = GM;
    end
end

if (~isempty(findobj('Tag','BundleFunc')))
    camUpVector = get(gca, 'CameraUpVector');
    camPosition = get(gca, 'CameraPosition');
    camTarget = get(gca, 'CameraTarget');
    camViewAngle = get(gca, 'CameraViewAngle');
    figurePosition = get(gcf, 'Position');
else
    figurePosition = [10, 10, 800, 800];
end

figure('Unit', 'pixel', 'Position', figurePosition);
set(gcf, 'ToolBar', 'none');
h = zeros(size(mesh_list));

for i=1:GroupSize
    color_data = repmat([0.9, 0.9, 0.8], mesh_list{i}.nV, 1);
    h(i) = subplot(DisplayLayout(1), DisplayLayout(2), DisplayOrder(i));
    mesh_list{i}.draw(struct('FaceColor', 'interp', 'FaceVertexCData', color_data, 'CDataMapping','scaled', 'EdgeColor', 'none', 'FaceAlpha',1,'AmbientStrength',0.3,'SpecularStrength',0.0));
    hold on;
    colormap jet(256);
    camlight('headlight');
    camlight(180,0);
    if strcmpi(options.names,'on')
        title(mesh_list{i}.Aux.name);
    end
    [IndsOnSource, ~] = GetLandmarks(mesh_list{i}.Aux.name,options.LandmarksPath,[options.MeshesPath mesh_list{i}.Aux.name options.MeshSuffix],options);
    draw_landmarks(mesh_list{i}.V, IndsOnSource);
end

if strcmpi(linkCamera, 'on')
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
end

if (exist('camUpVector', 'var'))
    set(gca, 'CameraUpVector', camUpVector);
    set(gca, 'CameraPosition', camPosition);
    set(gca, 'CameraTarget', camTarget);
    set(gca, 'CameraViewAngle', camViewAngle);
    set(gcf, 'Tag', 'BundleFunc');
end

end

function draw_landmarks(V,IndsOnSource,Type)

if nargin<3
    Type = 'Inds';
end

if strcmpi(Type,'Coords')
    Landmarks = IndsOnSource;
else
    if (size(V,1)==3)
        V = V';
    end
    Landmarks = V(IndsOnSource,:)';
end

NumLandmarks = max(size(Landmarks));
    
colmap =  [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
colmap = [colmap;colmap*0.5;colmap*0.25];

R = 0.015; % 0.025
CORR_draw_spheres(Landmarks',R,colmap(1:NumLandmarks,:));

end

