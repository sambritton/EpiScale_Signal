function signal_output = main_signaling(   arr_data_cell_index, arr_data_cell_pos_x, arr_data_cell_pos_y, ...
			arr_data_node_index, arr_data_node_pos_x, arr_data_node_pos_y)

% global thres_lateral thres_edge thres_nodes thres_neighbor
% thres_lateral = 0.64; % threshold for determining neighboring cells
% thres_edge = 0.005;
% thres_nodes = 0.75;  % need to adjust carefully to avoid triangles at junctions
% thres_neighbor = 0.5;
% 
% % load cell configuration from mechanical submodel
% %epi_nodes = load(['ExportCellProp_' num2str(index_couple) '.txt']);
% %epi_nodes(end,:) = [];
% 
% total_cell = max(cell_index);%index in c code has + 1 %max(epi_nodes(:,1))+1; 
% %cell_nodes = epi_nodes(total_cell+1:end,:);
% %cell_nodes = [cell_nodes(:,2:3), cell_nodes(:,1)];
% 
% cell_nodes = [cell_x, cell_y];%[node_pos_x, node_pos_y];
% for i = 1:total_cell
%     eval(['centroid_' num2str(i-1) '=[cell_center_x(i),cell_center_y(i)];']);
% end
% 
% dist_nodes = pdist2(cell_nodes(:,1:2),cell_nodes(:,1:2));
% % indicator of detection of common edge between neighboring cells
% common_edge_id = zeros(total_cell,total_cell);
% 
% for i = 0:total_cell-1
%     row_i = find(cell_nodes(:,3)==i);
%     dist_i = dist_nodes(row_i,:);
%     [row,col] = find(dist_i<thres_lateral);
%     neighbor_i = unique(cell_nodes(col,3)); %determine neighboring cells
%     neighbor_i = neighbor_i(neighbor_i~=i);
%     %construct common edge and indicator for cell-cell contact
%     cella = cell_nodes(cell_nodes(:,3)==i,1:2);
%     for j = 1:length(neighbor_i)
%         if common_edge_id(i+1,neighbor_i(j)+1)==0
%             cellb = cell_nodes(cell_nodes(:,3)==neighbor_i(j),1:2);
%             eval(['[intf_' num2str(i) '_' num2str(neighbor_i(j)) ', index_' num2str(i) '_' num2str(neighbor_i(j)) ', index_' num2str(neighbor_i(j)) '_' num2str(i) '] = detect_common_edge(cella,cellb);']);
%             common_edge_id(i+1,neighbor_i(j)+1)=1;
%             common_edge_id(neighbor_i(j)+1,i+1)=1;
%         end
%     end
% end
% 
% % detect vertices on each cell one by one
% for i = 0:total_cell-1
%     id = find(cell_nodes(:,3)==i);
%     newcell = cell_nodes(id,1:2);
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     delete_edge = [];
%     add_edge = [];
%     for j = 1:length(neighbor_i)
%         eval(['delete_edge = [delete_edge; index_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
%         if i<neighbor_i(j)
%             eval(['add_edge = [add_edge; intf_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
%         else
%             eval(['add_edge = [add_edge; intf_' num2str(neighbor_i(j)) '_' num2str(i) '];']);
%         end
%     end 
%     newcell(delete_edge,:) = [];
%     newcell = [newcell; add_edge];
%             
%     % sort the nodes counter-clockwisely
%     eval(['[newcell] = sort_counterclock(newcell,centroid_' num2str(i) ');']);
% %     plot(newcell(:,1),newcell(:,2),'*'); pause
% 
%     eval(['[vt_' num2str(i) '] = detect_vertices(newcell);']);
% end
% 
% % plot(vt_1(:,1),vt_1(:,2),'*','Color',[0 0 0]); pause
% 
% %choose the union set of vertices on common edge
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     for j = 1:length(neighbor_i)
%         if i<neighbor_i(j)
%             eval(['[C_keep,ia_keep,ib_keep] = intersect(vt_' num2str(neighbor_i(j)) ',intf_' num2str(i) '_' num2str(neighbor_i(j)) ',''rows'');']);
%         else
%             eval(['[C_keep,ia_keep,ib_keep] = intersect(vt_' num2str(neighbor_i(j)) ',intf_' num2str(neighbor_i(j)) '_' num2str(i) ',''rows'');']);
%         end
%         eval(['vt_' num2str(i) '= [vt_' num2str(i) ';C_keep];']);
%     end
%     eval(['[vt_' num2str(i) '] = unique(vt_' num2str(i) ',''rows'');']);
%     eval(['[vt_' num2str(i) '] = sort_tsp(vt_' num2str(i) ');']);
%     
% end
% 
% % for i = 0:total_cell-1
% %     eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''o'');']); 
% % end
% 
% 
% %replace points too close by the middle point on common edge
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     for j = 1:length(neighbor_i)
%         if neighbor_i(j)>i
%             eval(['[vt_contact,ia,ib] = intersect(vt_' num2str(i) ',vt_' num2str(neighbor_i(j)) ',''rows'');']);
%             eval(['[vt_contact] = sort_counterclock(vt_contact,centroid_' num2str(i) ');']); 
%         
%             eval(['[vt_contact_' num2str(i) '_' num2str(neighbor_i(j)) '] = combine_nodes(vt_contact);']);
%              
%             eval(['vt_' num2str(i) '(ia,:) = [];']); 
%             eval(['vt_' num2str(i) ' = [vt_' num2str(i) ';vt_contact_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
%             eval(['vt_' num2str(neighbor_i(j)) '(ib,:) = [];']);
%             eval(['vt_' num2str(neighbor_i(j)) '= [vt_' num2str(neighbor_i(j)) ';vt_contact_' num2str(i) '_' num2str(neighbor_i(j)) '];']);
%             eval(['[vt_' num2str(i) '] = sort_tsp(vt_' num2str(i) ');']);
%             eval(['[vt_' num2str(neighbor_i(j)) '] = sort_tsp(vt_' num2str(neighbor_i(j)) ');']);
%         end
%     end        
% end
% 
% % % visualize the triangular mesh
% % figure(1); hold on;
% % plot(cell_nodes(:,1),cell_nodes(:,2),'r.','MarkerSize',10);
% % for i = 0:total_cell-1
% %    eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''*'',''MarkerSize'',10);']);
% % end
% % 
% % for i = 58
% %    eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''.'',''MarkerSize'',20);']); 
% % end
% % 
% % for i = 0:total_cell-1
% %     eval(['temp = vt_' num2str(i) ';']);
% %    plot([temp(:,1); temp(1,1)],[temp(:,2); temp(1,2)],'b','LineWidth',2)
% %    for j = 1:size(temp,1)
% %        eval(['plot([centroid_' num2str(i) '(1) temp(' num2str(j) ',1)],[centroid_' num2str(i) '(2) temp(' num2str(j) ',2)],''b'',''LineWidth'',2);']);
% %    end
% % end
% 
% % romove vertices according to skewness
% for i = 0:total_cell-1
%     eval(['temp_vt = vt_' num2str(i) ';']);
%     eval(['temp_center = centroid_' num2str(i) ';']);
%     temp_skewness = zeros(size(temp_vt,1),1);
%     for ii = 1:size(temp_vt,1)
%         if ii==size(temp_vt,1)
%             temp_skewness(ii) = tri_skewness([temp_vt(ii,1); temp_vt(1,1); temp_center(1)],[temp_vt(ii,2); temp_vt(1,2); temp_center(2)]);
%         else
%             temp_skewness(ii) = tri_skewness([temp_vt(ii:ii+1,1); temp_center(1)],[temp_vt(ii:ii+1,2); temp_center(2)]);
%         end
%     end
%     [val,ind] = max(temp_skewness);
%     while val > 0.75
%         if ind>1
%             L = temp_skewness(ind-1);
%         else
%             L = temp_skewness(end);
%         end
%         if ind<length(temp_skewness)
%             R = temp_skewness(ind+1);
%         else
%             R = temp_skewness(1);
%         end
%         if L>=R
%             ind_skewness = ind;
%         else
%             ind_skewness = (ind+1)*(ind<length(temp_skewness))+1*(ind==length(temp_skewness));
%         end
%         eval(['pt_skewness = vt_' num2str(i) '(ind_skewness,:);']);
%         eval(['vt_' num2str(i) '(ind_skewness,:)=[];']);
%         neighbor_i = find(common_edge_id(i+1,:));
%         neighbor_i = neighbor_i - 1;
%         for j = 1:length(neighbor_i)
%             eval(['[id1,ind1] = ismember(pt_skewness,vt_' num2str(neighbor_i(j)) ',''rows'');']);
%             if id1 == 1
%                 eval(['vt_' num2str(neighbor_i(j)) '(ind1,:) = [];']);
%             end
%         end
%         eval(['temp_vt = vt_' num2str(i) ';']);
%         temp_skewness = zeros(size(temp_vt,1),1);
%         for ii = 1:size(temp_vt,1)
%             if ii==size(temp_vt,1)
%                 temp_skewness(ii) = tri_skewness([temp_vt(ii,1); temp_vt(1,1); temp_center(1)],[temp_vt(ii,2); temp_vt(1,2); temp_center(2)]);
%             else
%                 temp_skewness(ii) = tri_skewness([temp_vt(ii:ii+1,1); temp_center(1)],[temp_vt(ii:ii+1,2); temp_center(2)]);
%             end
%         end
%         [val,ind] = max(temp_skewness);
%     end
% end
% 
%        
% % % visualize the triangular mesh
% % figure(2); hold on;
% % plot(cell_nodes(:,1),cell_nodes(:,2),'r.','MarkerSize',10);
% % for i = 0:total_cell-1
% %    eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''*'',''MarkerSize'',10);']); 
% % end
% % 
% % for i = 0:total_cell-1
% %     eval(['temp = vt_' num2str(i) ';']);
% %    plot([temp(:,1); temp(1,1)],[temp(:,2); temp(1,2)],'b','LineWidth',2)
% %    for j = 1:size(temp,1)
% %        eval(['plot([centroid_' num2str(i) '(1) temp(' num2str(j) ',1)],[centroid_' num2str(i) '(2) temp(' num2str(j) ',2)],''b'',''LineWidth'',2);']);
% %    end
% % end
% 
% 
% % replace junction points among multiple(>=3) contacting cells by centroids
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     for j = 1:length(neighbor_i)
%         neighbor_j = find(common_edge_id(neighbor_i(j)+1,:));
%         neighbor_j = neighbor_j - 1;
%         neighbor_ij = intersect(neighbor_i,neighbor_j);
%         for k = 1:length(neighbor_ij)
%             c1 = i;
%             c2 = neighbor_i(j);
%             c3 = neighbor_ij(k);
%             if c1<c2 && c2<c3
%                 eval(['e12 = intersect(vt_' num2str(c1) ',vt_' num2str(c2) ',''rows'');']);
%                 eval(['e13 = intersect(vt_' num2str(c1) ',vt_' num2str(c3) ',''rows'');']);
%                 eval(['e23 = intersect(vt_' num2str(c2) ',vt_' num2str(c3) ',''rows'');']);
%                 if ~isempty(e12) && ~isempty(e13) && ~isempty(e23)
%                     [v12,v13,v23] = find_inner_tri(e12,e13,e23);
%                     if ~isempty(v12) && ~isempty(v13) && ~isempty(v23)
%                         eval(['[id1_12,ind1_12] = ismember(v12,vt_' num2str(c1) ',''rows'');']);
%                         eval(['[id2_12,ind2_12] = ismember(v12,vt_' num2str(c2) ',''rows'');']);
%                         eval(['[id1_13,ind1_13] = ismember(v13,vt_' num2str(c1) ',''rows'');']);
%                         eval(['[id3_13,ind3_13] = ismember(v13,vt_' num2str(c3) ',''rows'');']);
%                         eval(['[id2_23,ind2_23] = ismember(v23,vt_' num2str(c2) ',''rows'');']);
%                         eval(['[id3_23,ind3_23] = ismember(v23,vt_' num2str(c3) ',''rows'');']);
%                         centroid123 = (v12+v13+v23)/3;
%                         
%                         eval(['vt_' num2str(c1) '(ind1_12,:) = centroid123;']);
%                         eval(['vt_' num2str(c2) '(ind2_12,:) = centroid123;']);
%                         eval(['vt_' num2str(c1) '(ind1_13,:) = centroid123;']);
%                         eval(['vt_' num2str(c3) '(ind3_13,:) = centroid123;']);
%                         eval(['vt_' num2str(c2) '(ind2_23,:) = centroid123;']);
%                         eval(['vt_' num2str(c3) '(ind3_23,:) = centroid123;']);
%                         eval(['[vt_' num2str(c1) '] = unique(vt_' num2str(c1) ',''rows'');']);
%                         eval(['[vt_' num2str(c2) '] = unique(vt_' num2str(c2) ',''rows'');']);
%                         eval(['[vt_' num2str(c3) '] = unique(vt_' num2str(c3) ',''rows'');']);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % check whether vertices that are not the boundary of the tissue are shared
% % by more than one cell; if not, add the vertice to the neighboring cells
% % with distance below some threshold
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     vet_neighbor = [];
%     for j = 1:length(neighbor_i)
%         eval(['vet_neighbor = [vet_neighbor; vt_' num2str(neighbor_i(j)) '];']);
%     end
%     eval(['num_vet_i = size(vt_' num2str(i) ',1);']);
%     for j = 1:num_vet_i
%        eval(['temp = vt_' num2str(i) '(j,:); ']);
%        if ~ismember(temp,vet_neighbor)
%            for k = 1:length(neighbor_i)
%                id = find(cell_nodes(:,3)==neighbor_i(k));
%                neighborcell = cell_nodes(id,1:2);
%                dist = pdist2(temp,neighborcell);
%                if min(dist) < thres_neighbor
%                    eval(['vt_' num2str(neighbor_i(k)) ' = [vt_' num2str(neighbor_i(k)) '; temp];']);
%                end
%            end
%        elseif sum(vet_neighbor == temp)==1
%            for k = 1:length(neighbor_i)
%                id = find(cell_nodes(:,3)==neighbor_i(k));
%                neighborcell = cell_nodes(id,1:2);
%                dist = pdist2(temp,neighborcell);
%                if min(dist) < thres_neighbor
%                    eval(['vt_' num2str(neighbor_i(k)) ' = [vt_' num2str(neighbor_i(k)) '; temp];']);
%                end
%            end
%        end
%     end
% end
% % for i = 0:total_cell-1
% %     eval(['[vt_' num2str(i) '] = unique(vt_' num2str(i) ',''rows'');']);
% %     vt_num = eval(['size(vt_' num2str(i) ',1)']);
% %     eval(['vt_dist = sqrt((vt_' num2str(i) '(:,1)-[vt_' num2str(i) '(2:end,1); vt_' num2str(i) '(1,1)]).^2+(vt_' num2str(i) '(:,2)-[vt_' num2str(i) '(2:end,2); vt_' num2str(i) '(1,2)]).^2);']);
% %     [val,ind] = min(vt_dist);
% %     while val<thres_nodes
% %         if ind < vt_num
% %             eval(['pt1 = vt_' num2str(i) '(ind,:);']);
% %             eval(['pt2 = vt_' num2str(i) '(ind+1,:);']);
% %             eval(['mid_pt = ( vt_' num2str(i) '(ind,:)+vt_' num2str(i) '(ind+1,:) )/2;']);
% %             eval(['vt_' num2str(i) '(ind,:) = mid_pt;']);
% %             eval(['vt_' num2str(i) '(ind+1,:) = [];']);
% %         else
% %             eval(['pt1 = vt_' num2str(i) '(end,:);']);
% %             eval(['pt2 = vt_' num2str(i) '(1,:);']);
% %             eval(['mid_pt = ( vt_' num2str(i) '(end,:)+vt_' num2str(i) '(1,:) )/2;']);
% %             eval(['vt_' num2str(i) '(end,:) = mid_pt;']);
% %             eval(['vt_' num2str(i) '(1,:) = [];']);
% %         end
% %         neighbor_i = find(common_edge_id(i+1,:));
% %         neighbor_i = neighbor_i - 1;
% %         for j = 1:length(neighbor_i)
% %             eval(['[id1,ind1] = ismember(pt1,vt_' num2str(neighbor_i(j)) ',''rows'');']);
% %             if id1 == 1
% %                 eval(['vt_' num2str(neighbor_i(j)) '(ind1,:) = mid_pt;']);
% %             end
% %             eval(['[id2,ind2] = ismember(pt2,vt_' num2str(neighbor_i(j)) ',''rows'');']);
% %             if id2 == 1
% %                 eval(['vt_' num2str(neighbor_i(j)) '(ind2,:) = mid_pt;']);
% %             end
% %         end
% %         vt_num = vt_num -1;
% %         eval(['vt_dist = sqrt((vt_' num2str(i) '(:,1)-[vt_' num2str(i) '(2:end,1); vt_' num2str(i) '(1,1)]).^2+(vt_' num2str(i) '(:,2)-[vt_' num2str(i) '(2:end,2); vt_' num2str(i) '(1,2)]).^2);']);
% %         [val,ind] = min(vt_dist);
% %     end
% % end
% 
% % avoid repeated points and set the orientation of the vertices
% for i = 0:total_cell-1
%     eval(['[vt_' num2str(i) '] = unique(vt_' num2str(i) ',''rows'');']);
%     eval(['[vt_' num2str(i) '] = sort_tsp(vt_' num2str(i) ');']);
% end
% 
% % romove vertices according to skewness
% for i = 0:total_cell-1
%     eval(['temp_vt = vt_' num2str(i) ';']);
%     eval(['temp_center = centroid_' num2str(i) ';']);
%     temp_skewness = zeros(size(temp_vt,1),1);
%     for ii = 1:size(temp_vt,1)
%         if ii==size(temp_vt,1)
%             temp_skewness(ii) = tri_skewness([temp_vt(ii,1); temp_vt(1,1); temp_center(1)],[temp_vt(ii,2); temp_vt(1,2); temp_center(2)]);
%         else
%             temp_skewness(ii) = tri_skewness([temp_vt(ii:ii+1,1); temp_center(1)],[temp_vt(ii:ii+1,2); temp_center(2)]);
%         end
%     end
%     [val,ind] = max(temp_skewness);
%     while val > 0.75
%         if ind>1
%             L = temp_skewness(ind-1);
%         else
%             L = temp_skewness(end);
%         end
%         if ind<length(temp_skewness)
%             R = temp_skewness(ind+1);
%         else
%             R = temp_skewness(1);
%         end
%         if L>=R
%             ind_skewness = ind;
%         else
%             ind_skewness = (ind+1)*(ind<length(temp_skewness))+1*(ind==length(temp_skewness));
%         end
%         eval(['pt_skewness = vt_' num2str(i) '(ind_skewness,:);']);
%         eval(['vt_' num2str(i) '(ind_skewness,:)=[];']);
%         neighbor_i = find(common_edge_id(i+1,:));
%         neighbor_i = neighbor_i - 1;
%         for j = 1:length(neighbor_i)
%             eval(['[id1,ind1] = ismember(pt_skewness,vt_' num2str(neighbor_i(j)) ',''rows'');']);
%             if id1 == 1
%                 eval(['vt_' num2str(neighbor_i(j)) '(ind1,:) = [];']);
%             end
%         end
%         eval(['temp_vt = vt_' num2str(i) ';']);
%         temp_skewness = zeros(size(temp_vt,1),1);
%         for ii = 1:size(temp_vt,1)
%             if ii==size(temp_vt,1)
%                 temp_skewness(ii) = tri_skewness([temp_vt(ii,1); temp_vt(1,1); temp_center(1)],[temp_vt(ii,2); temp_vt(1,2); temp_center(2)]);
%             else
%                 temp_skewness(ii) = tri_skewness([temp_vt(ii:ii+1,1); temp_center(1)],[temp_vt(ii:ii+1,2); temp_center(2)]);
%             end
%         end
%         [val,ind] = max(temp_skewness);
%     end
% end
% 
% % add junction points among 3 contacting cells by centroids
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     for j = 1:length(neighbor_i)
%         neighbor_j = find(common_edge_id(neighbor_i(j)+1,:));
%         neighbor_j = neighbor_j - 1;
%         neighbor_ij = intersect(neighbor_i,neighbor_j);
%         for k = 1:length(neighbor_ij)
%             c1 = i;
%             c2 = neighbor_i(j);
%             c3 = neighbor_ij(k);
%             if c1<c2 && c2<c3
%                 eval(['e12 = intersect(vt_' num2str(c1) ',vt_' num2str(c2) ',''rows'');']);
%                 eval(['e13 = intersect(vt_' num2str(c1) ',vt_' num2str(c3) ',''rows'');']);
%                 eval(['e23 = intersect(vt_' num2str(c2) ',vt_' num2str(c3) ',''rows'');']);
%                 if ~isempty(e12) && ~isempty(e13) && ~isempty(e23)
%                     [v12,v13,v23] = find_inner_tri(e12,e13,e23);
%                     if ~isempty(v12) && ~isempty(v13) && ~isempty(v23)
%                         eval(['[id1_12,ind1_12] = ismember(v12,vt_' num2str(c1) ',''rows'');']);
%                         eval(['[id2_12,ind2_12] = ismember(v12,vt_' num2str(c2) ',''rows'');']);
%                         eval(['[id1_13,ind1_13] = ismember(v13,vt_' num2str(c1) ',''rows'');']);
%                         eval(['[id3_13,ind3_13] = ismember(v13,vt_' num2str(c3) ',''rows'');']);
%                         eval(['[id2_23,ind2_23] = ismember(v23,vt_' num2str(c2) ',''rows'');']);
%                         eval(['[id3_23,ind3_23] = ismember(v23,vt_' num2str(c3) ',''rows'');']);
%                         centroid123 = (v12+v13+v23)/3;
%                         
%                         eval(['vt_' num2str(c1) ' = [vt_' num2str(c1) '; centroid123];']);
%                         eval(['vt_' num2str(c2) ' = [vt_' num2str(c2) '; centroid123];']);
%                         eval(['vt_' num2str(c3) ' = [vt_' num2str(c3) '; centroid123];']);
%                     end
%                 end
%             end
%         end
%     end
% end
% 
% % add junction points among 4 contacting cells by centroids
% for i = 0:total_cell-1
%     neighbor_i = find(common_edge_id(i+1,:));
%     neighbor_i = neighbor_i - 1;
%     for j = 1:length(neighbor_i)
%         neighbor_j = find(common_edge_id(neighbor_i(j)+1,:));
%         neighbor_j = neighbor_j - 1;
%         for k = 1:length(neighbor_j)
%             neighbor_k = find(common_edge_id(neighbor_j(k)+1,:));
%             neighbor_k = neighbor_k - 1;
%             neighbor_ki = intersect(neighbor_i,neighbor_k);
%             for l = 1:length(neighbor_ki)
%                 c1 = i;
%                 c2 = neighbor_i(j);
%                 c3 = neighbor_j(k);
%                 c4 = neighbor_ki(l);
%                 
%                 eval(['e12 = intersect(vt_' num2str(c1) ',vt_' num2str(c2) ',''rows'');']);
%                 eval(['e23 = intersect(vt_' num2str(c2) ',vt_' num2str(c3) ',''rows'');']);
%                 eval(['e34 = intersect(vt_' num2str(c3) ',vt_' num2str(c4) ',''rows'');']);
%                 eval(['e41 = intersect(vt_' num2str(c4) ',vt_' num2str(c1) ',''rows'');']);
%                 
%                 if ~isempty(e12) && ~isempty(e23) && ~isempty(e34) && ~isempty(e41)
%                     common_123 = intersect(e12,e23,'rows');
%                     common_234 = intersect(e23,e34,'rows');
%                     common_341 = intersect(e34,e41,'rows');
%                     common_412 = intersect(e41,e12,'rows');
%                     junction_1234 = [common_123; common_234; common_341; common_412];
%                     if isempty(junction_1234)
%                         [v12,v23,v34,v41] = find_inner_quad(e12,e23,e34,e41);
%                         centroid1234 = (v12+v23+v34+v41)/4;
%                         if min(pdist2(centroid1234,epi_nodes(1:total_cell,2:3)))>1
%                             eval(['vt_' num2str(c1) ' = [vt_' num2str(c1) '; centroid1234];']);
%                             eval(['vt_' num2str(c2) ' = [vt_' num2str(c2) '; centroid1234];']);
%                             eval(['vt_' num2str(c3) ' = [vt_' num2str(c3) '; centroid1234];']);
%                             eval(['vt_' num2str(c4) ' = [vt_' num2str(c4) '; centroid1234];']);
%                         end
%                     end
%                 end
%                 
%             end
%         end
%     end
% end
% 
% % avoid repeated points and set the orientation of the vertices
% for i = 0:total_cell-1
%     eval(['[vt_' num2str(i) '] = unique(vt_' num2str(i) ',''rows'');']);
%     eval(['[vt_' num2str(i) '] = sort_tsp(vt_' num2str(i) ');']);
% end
% 
% 
% % figure(3); hold on;
% % plot(cell_nodes(:,1),cell_nodes(:,2),'r.','MarkerSize',10);
% % for i = 0:total_cell-1
% %    eval(['plot(vt_' num2str(i) '(:,1),vt_' num2str(i) '(:,2),''*'',''MarkerSize'',10);']); 
% % end
% % for i = 0:total_cell-1
% %     eval(['temp = vt_' num2str(i) ';']);
% %    plot([temp(:,1); temp(1,1)],[temp(:,2); temp(1,2)],'b','LineWidth',2)
% %    for j = 1:size(temp,1)
% %        eval(['plot([centroid_' num2str(i) '(1) temp(' num2str(j) ',1)],[centroid_' num2str(i) '(2) temp(' num2str(j) ',2)],''b'',''LineWidth'',2);']);
% %    end
% % end
% 
% 
% %saveas(figure(1),'triangle_mesh.fig'); 
% 
% 
% % solve u_t = D(u_xx+u_yy)-du+source on the triangular mesh by
% % approximating diffusion by passive transport and imposing source function
% % at cells close to 25
% mesh_size = 0;
% vt_all = [];
% for i = 0:total_cell-1
%     eval(['num_tri_' num2str(i) ' = size(vt_' num2str(i) ',1);']);
%     eval(['vt_all = [vt_all; vt_' num2str(i) '];'])
%     eval(['mesh_size = mesh_size + num_tri_' num2str(i) ';']);
% end
% u = zeros(mesh_size,1);
% % construct a matrix to denote neighboring triangles
% NI_mat = [];
% % construct a matrix to store contact length
% A_mat = [];
% % intracellular neighbors
% for i = 0:total_cell-1
%     eval(['temp = ones(num_tri_' num2str(i) ',1);']);
%     eval(['NI_i = spdiags(temp,1,num_tri_' num2str(i) ',num_tri_' num2str(i) ')+spdiags(temp,-1,num_tri_' num2str(i) ',num_tri_' num2str(i) ');']);
%     eval(['NI_i(1,num_tri_' num2str(i) ') = 1;']);
%     eval(['NI_i(num_tri_' num2str(i) ',1) = 1;']);
%     NI_mat = blkdiag(NI_mat,NI_i);
%     
%     eval(['edge_vec = pdist2(centroid_' num2str(i) ',vt_' num2str(i) ');']);
%     edge_vec = [edge_vec(2:end), edge_vec(1)];
%     edge_vec = edge_vec';
%     eval(['A_' num2str(i) ' = spdiags(edge_vec,1,num_tri_' num2str(i) ',num_tri_' num2str(i) ') + spdiags(edge_vec,-1,num_tri_' num2str(i) ',num_tri_' num2str(i) ');']);
%     eval(['A_' num2str(i) '(1,end) = edge_vec(end);']);
%     eval(['A_' num2str(i) '(end,1) = edge_vec(end);']);
%     eval(['A_mat = blkdiag(A_mat,A_' num2str(i) ');']); 
% end
% % intercellular neighbors
% % be aware of the junction points: ind1 and ind2 may have length greater
% % than 2
% for i = 0:total_cell-1
%     eval(['num_tri = num_tri_' num2str(i) ';']);
%     for j = 1:num_tri
%         if j == num_tri
%             jnext = 1;
%         else
%             jnext = j+1;
%         end
%         eval(['[ind1] = find(sum(abs(vt_all-repmat(vt_' num2str(i) '(j,:),mesh_size,1)),2)==0);']);
%         eval(['[ind2] = find(sum(abs(vt_all-repmat(vt_' num2str(i) '(jnext,:),mesh_size,1)),2)==0);']);
%         if length(ind1)>1 && length(ind2)>1
%             eval(['Ledge = pdist(vt_' num2str(i) '([j jnext],:));']);
%             % find two pairs with minimum differences corresponding to
%             % contacting edges
%             diff = repmat(ind1,1,length(ind2))-repmat(ind2',length(ind1),1);
%             diff = abs(diff);
%             [val,I2] = min(diff,[],2);
%             [val,I1] = sort(val);
%             common_edge = [ind1(I1(1)),ind2(I2(I1(1)));ind1(I1(2)),ind2(I2(I1(2)))];
%             if abs(common_edge(1,1)-common_edge(1,2)) ==1
%                 tri_1 = min(common_edge(1,1),common_edge(1,2));
%             else
%                 tri_1 = max(common_edge(1,1),common_edge(1,2));
%             end
%             if abs(common_edge(2,1)-common_edge(2,2)) ==1
%                 tri_2 = min(common_edge(2,1),common_edge(2,2));
%             else
%                 tri_2 = max(common_edge(2,1),common_edge(2,2));
%             end
%             NI_mat(tri_1,tri_2) = 1;
%             NI_mat(tri_2,tri_1) = 1;
%             A_mat(tri_1,tri_2) = Ledge;
%             A_mat(tri_2,tri_1) = Ledge;
%         end
%     end
% end
% 
% % centroids of each triangle
% mesh_centroids = zeros(mesh_size,2);
% num_tri_pre = 0;
% for i = 0:total_cell-1
%     eval(['num_tri = num_tri_' num2str(i) ';']);
%     for j = 1:num_tri
%         if j == num_tri
%             jnext = 1;
%         else
%             jnext = j+1;
%         end
%         eval(['mesh_centroids(num_tri_pre+j,:) =  sum([vt_' num2str(i) '([j jnext],:);centroid_' num2str(i) '])/3;']);
%     end
%     num_tri_pre = num_tri_pre + num_tri; 
% end
% L_mat = pdist2(mesh_centroids,mesh_centroids);
% L_mat = L_mat + 10000*eye(mesh_size); %avoid singularities
% 
% % source cells: located within 12% of the total tissue size around the
% % midline
% tissue_centroid = zeros(1,2);
% for i = 1:total_cell
%     eval(['tissue_centroid = tissue_centroid + centroid_' num2str(i-1) ';']);
% end
% tissue_centroid = tissue_centroid/total_cell;
% tissue_r = sqrt( sum((tissue_centroid-centroid_0).^2) );
% for i = 2:total_cell
%     eval(['temp = sqrt( sum((tissue_centroid-centroid_' num2str(i-1) ').^2) );']);
%     if temp > tissue_r
%         tissue_r = temp;
%     end
% end
%     
% [ind] = find(abs(cell_nodes(:,1)-tissue_centroid(1))<0.12*tissue_r);
% source_cellid = unique(cell_nodes(ind,3));
% source_cellid = sort(source_cellid,'ascend');
% u_source = [];
% for i = 0:total_cell-1
%     temp = ismember(i,source_cellid);
%     if temp == 1
%         eval(['u_source = [u_source; ones(num_tri_' num2str(i) ',1)];']);
%     else
%         eval(['u_source = [u_source; zeros(num_tri_' num2str(i) ',1)];']);
%     end
% end
% 
% Dpp_mat = zeros(mesh_size,4);
% dt = 0.002;
% 
% % keyboard
% 
% for iter = 1:100
%     [frhs] = Dpp_signaling(Dpp_mat,A_mat,L_mat,NI_mat,u_source);
%     temp = Dpp_mat + dt * frhs;
%     
%     [frhs] = Dpp_signaling(temp,A_mat,L_mat,NI_mat,u_source);
%     Dpp_mat = Dpp_mat/2 + temp/2 + dt/2 * frhs; 
% %     display(iter)
% end
% 
% Dpp = Dpp_mat(:,1);
% Dpp_cell = zeros(total_cell,1);
% Tkv = Dpp_mat(:,2);
% Tkv_cell = zeros(total_cell,1);
% DT = Dpp_mat(:,3);
% DT_cell = zeros(total_cell,1);
% pMad = Dpp_mat(:,4);
% pMad_cell = zeros(total_cell,1);
% % figure(2); 
% % % subplot(2,2,1); 
% % hold on;
% num_tri_pre=0;
% for i = 0:total_cell-1
%     eval(['num_tri = num_tri_' num2str(i) ';']);
%     for j = 1:num_tri
%         if j == num_tri
%             jnext = 1;
%         else
%              jnext = j+1;
%          end
% %         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],Dpp(num_tri_pre+j),''LineStyle'',''none'');']);
%      end
%      Dpp_cell(i+1) = sum(Dpp(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
%      num_tri_pre = num_tri_pre + num_tri;
%  end
% % colorbar;
% % subplot(2,2,2); hold on;
%  num_tri_pre=0;
%  for i = 0:total_cell-1
%      eval(['num_tri = num_tri_' num2str(i) ';']);
%      for j = 1:num_tri
%          if j == num_tri
%              jnext = 1;
%          else
%              jnext = j+1;
%          end
% %         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],Tkv(num_tri_pre+j));']);
%      end
%      Tkv_cell(i+1) = sum(Tkv(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
%      num_tri_pre = num_tri_pre + num_tri;
%  end
% % colorbar;
% % subplot(2,2,3); hold on;
%  num_tri_pre=0;
%  for i = 0:total_cell-1
%      eval(['num_tri = num_tri_' num2str(i) ';']);
%      for j = 1:num_tri
%          if j == num_tri
%              jnext = 1;
%          else
%              jnext = j+1;
%          end
% %         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],DT(num_tri_pre+j));']);
%      end
%      DT_cell(i+1) = sum(DT(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
%      num_tri_pre = num_tri_pre + num_tri;
%  end
% % colorbar;
% % subplot(2,2,4); hold on;
%  num_tri_pre=0;
%  for i = 0:total_cell-1
%      eval(['num_tri = num_tri_' num2str(i) ';']);
%      for j = 1:num_tri
%          if j == num_tri
%              jnext = 1;
%          else
%              jnext = j+1;
%          end
% %         eval(['fill([vt_' num2str(i) '(j,1) vt_' num2str(i) '(jnext,1) centroid_' num2str(i) '(1)],[vt_' num2str(i) '(j,2) vt_' num2str(i) '(jnext,2) centroid_' num2str(i) '(2)],pMad(num_tri_pre+j));']);
%      end
%      pMad_cell(i+1) = sum(pMad(num_tri_pre+1:num_tri_pre+num_tri))/num_tri;
%      num_tri_pre = num_tri_pre + num_tri;
%  end
% % colorbar;
% % saveas(figure(2),'Dpp_signaling.fig');
% % figure(2);
% % plot(cell_nodes(:,1),cell_nodes(:,2),'.','Color',[1 1 1],'MarkerSize',12);
% % ch = colorbar('Ticks',[0.01 0.07],'TickLabels',{'low','high'});
% % set(ch,'FontSize',25);
% % ylabel(ch,'Dpp signaling')
% % filename = ['Dpp_cell_T' num2str(index_couple) '.txt'];
% % %fileID = fopen(filename,'w');
% % fileID = fopen('tmp.txt','w+');
% % fprintf(fileID,'%.4f\n',Dpp_cell);
% % fclose(fileID);
% % movefile ('tmp.txt',filename) ; 
Dpp_cell = [1,2,3];

signal_output = Dpp_cell;

