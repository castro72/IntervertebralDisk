function valid = checkUniqueDim(Xg)
valid = true;
for i = 1:size(Xg,2)
    Xr = Xg(:,i);
    if numel(unique(Xr))==1
        valid = false;
        break;
    end
end