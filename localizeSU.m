% localize the neuron

cid_list=h5read('transient_6.hdf5','/cluster_id');
path_list=h5read('transient_6.hdf5','/path');
%reg_list=h5read('transient_6.hdf5','/reg'); % positive control

pathFullList = h5read('path2tid.hdf5','/path');
midList = h5read('path2tid.hdf5','/mid');
tidList = h5read('path2tid.hdf5','/tid');
depthList = h5read('path2tid.hdf5','/depth');

% in case all tracks were labeled imec0
path_list(cid_list>10000) = cellfun(@(x) strrep(x,'imec0','imec1'),path_list(cid_list>10000),'UniformOutput',false);
cid_list(cid_list>10000) = cid_list(cid_list>10000) - 10000;


supool=1:length(cid_list);

done=[];
error_list=cell(0);
coord = nan(length(cid_list),3);
% directory of reference atlas files
annotation_volume_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\structure_tree_safe_2017.csv';
% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

prevMid = 0;
for i=1:length(supool)
    if ismember(supool(i),done)
        continue
    end
    folder=deblank(path_list{supool(i)});
    currMid = midList(contains(pathFullList,folder));
    currTid = tidList(contains(pathFullList,folder));
    depth = depthList(contains(pathFullList,folder));
    if ~isempty(currMid) && isfile(fullfile('H:\NP histology\probe_points',sprintf('%d.mat',currMid))) && currTid>0 
        if currMid ~= prevMid
            load(fullfile('H:\NP histology\probe_points',sprintf('%d.mat',currMid)));
            prevMid = currMid;
        end
        cIds = cid_list(startsWith(path_list,folder));
        fid = fopen(fullfile(pathFullList{contains(pathFullList,folder)},'cluster_info.tsv'));
        T = textscan(fid,'%d %f %f %s %f %d %f %s %s %d %d', 'HeaderLines', 1); 
        fclose(fid);
        dist2TipList = T{7};
        cInfoId = T{1};
        [~,tf,cidx] = intersect(double(cIds),double(cInfoId),'stable');
        realDepth = double(depth) - dist2TipList(cidx);
        try
            curr_probePoints = pointList.pointList{currTid,1}(:, [3 2 1]);
        catch
            disp('track id exceed track number');
            error_list(end+1,:) = {folder,4};
            continue;
        end
        % get line of best fit through points
        % m is the mean value of each dimension; p is the eigenvector for largest eigenvalue
        [m,p,s] = best_fit_line(curr_probePoints(:,1), curr_probePoints(:,2), curr_probePoints(:,3));
        if isnan(m(1))
            disp('no points found current probe')
            error_list(end+1,:) = {folder,3};
            continue;
        end
        % ensure proper orientation: want 0 at the top of the brain and positive distance goes down into the brain
        if p(2)<0
            p = -p;
        end
        % determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
        ann = 10;
        out_of_brain = false;
        while ~(ann==1 && out_of_brain) % && distance_stepped > .5*active_probe_length)
            m = m-p; % step 10um, backwards up the track
            ann = av(round(m(1)),round(m(2)),round(m(3))); %until hitting the top
            if strcmp(st.safe_name(ann), 'root')
                % make sure this isn't just a 'root' area within the brain
                m_further_up = m - p*20; % is there more brain 200 microns up along the track?
                ann_further_up = av(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
                if strcmp(st.safe_name(ann_further_up), 'root')
                    out_of_brain = true;
                end
            end
        end
        dest = [m(1)+p(1)*realDepth/10,m(2)+p(2)*realDepth/10,m(3)+p(3)*realDepth/10];
        try
            trackLoc = coord(startsWith(path_list,folder),:);
            trackLoc(tf,:) = dest;
            coord(startsWith(path_list,folder),:) = trackLoc;
        catch
            disp('SU number do not match');
            error_list(end+1,:) = {folder,2};
            continue;
        end
        done = [done;find(startsWith(path_list,folder))];
    else
        disp('Undefined track');
        error_list(end+1,:) = {folder,1};
        continue;
    end
end
 save('E:\prJ\neuropixels\histology location analysis\sucoords318.mat','coord','-v7.3');       

 
%% test

% acryms = cell(length(cid_list),1);
% for i = 1:length(cid_list)
%     x = coord(i,1);
%     y = coord(i,2);
%     z = coord(i,3);
%     if isnan(x)
%         continue;
%     end
%     if x>0 && x<=size(av,1) &&...
%        y>0 && y<=size(av,2) &&...
%        z>0 && z<=size(av,3)
%         ann = av(ceil(x), ceil(y), ceil(z));
%         acryms{i} = st.acronym{ann};
%     end
% end
% t = [acryms,reg_list];

