function [vt] = detect_vertices(cell)
global thres_edge
thres_edge = polyarea(cell(:,1),cell(:,2))/320;
vt = cell(1,:);
istart = 1;
iend = 1;
A = 0;
while iend<=size(cell,1)-2
    while A < thres_edge && iend<=size(cell,1)-1
        iend = iend + 1;
        A = polyarea(cell(istart:iend,1),cell(istart:iend,2));
       
    end
%     keyboard
    vt = [vt; cell(iend-1,:)];
    istart = iend-1;
    iend = iend -1;
    A = 0;
end
vt(end,:)=[];