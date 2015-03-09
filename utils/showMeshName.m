function output_txt = showMeshName(obj, event_obj, Y, Names)
%SHOWMESHNAME Summary of this function goes here
%   Detailed explanation goes here

pos = get(event_obj,'Position');
res = sum((Y-repmat(pos, size(Y,1), 1)).^2,2);
output_txt = num2str(Names{res < 1e-10});

end

