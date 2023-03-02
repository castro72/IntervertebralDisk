function data_out = resize_tendonImg(data_in,scales_in,iso_res)

cols = size(data_in,2);
rows = size(data_in,1);
tpts = size(data_in,3);

xaxis = (0:(cols - 1)) .* scales_in(2);
yaxis = (0:(rows - 1)) .* scales_in(1);

if isempty(iso_res)
    new_res = min(scales_in(1:2));
else
    new_res = iso_res;
end

xaxNew = (0:new_res:xaxis(end));
yaxNew = (0:new_res:yaxis(end));

[Xorig,Yorig] = meshgrid(xaxis,yaxis);
[Xnew, Ynew ] = meshgrid(xaxNew,yaxNew);

data_out = zeros(numel(yaxNew),numel(xaxNew),tpts,'uint8');

fprintf('Loading Data:\n');
checkpoints = round([0.1:0.1:1.0] .* tpts);
checklabels = [10:10:100];
count = 1;
for t = 1:tpts
    data_out(:,:,t) = uint8(interp2(Xorig,Yorig,single(data_in(:,:,t)),Xnew,Ynew));
    if t == checkpoints(count)
        fprintf('%g %%\n',checklabels(count));
        count = count + 1;
    end
end
% close(h1);