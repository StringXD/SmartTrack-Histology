function [sortedMatrix,order] = sortByPeakTimeAsc(oriMatrix,orderMatrix,dim,ordertype)
% sort neurons by the occurence of peak activity
if strcmp(ordertype,'ascend')
    [~,I] = max(orderMatrix,[],dim);
elseif strcmp(ordertype,'descend')
    [~,I] = min(orderMatrix,[],dim);
end
[~,II] = sort(I,'ascend');
if dim == 1
    sortedMatrix = oriMatrix(:,II);
elseif dim == 2
    sortedMatrix = oriMatrix(II,:);
end
if nargout > 1
    order = II;
end