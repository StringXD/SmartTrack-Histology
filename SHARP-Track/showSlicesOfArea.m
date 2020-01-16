% show slices and points given the brain area
%**************************************************************************
function[] = showSlicesOfArea(string)
% string is a string including acronym of brain areas, / as separator

% location of object points
pointListFolder = 'H:\NP histology\probe_points';

% location of processed folders
processedFolderRoot = 'H:\NP histology\processed';

% location of trackTable
trackTable = 'C:\Users\xd\Desktop\NP tracks validated.csv';

% location of result folders
rtsFolderRoot = 'H:\NP histology\results';

% location used to place slices of certain area
saveLocation = 'H:\NP histology\slices points';

% directory of reference atlas files
annotation_volume_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\structure_tree_safe_2017.csv';
template_volume_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\template_volume_10um.npy';
% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
    tv = readNPY(template_volume_location);
end

% load the track table
T = readtable(trackTable);

% separate the string
targetAeras = regexp(string,'\w+','match');

for currAeraIdx = 1:length(targetAeras)
    cd(saveLocation);
    if ~exist(targetAeras{currAeraIdx},'dir')
        mkdir(targetAeras{currAeraIdx});
    end
    cd(targetAeras{currAeraIdx});
    isA = isAreaOrContains(st.acronym, targetAeras{currAeraIdx}, st);
    avList = st(isA,:).sphinx_id;
    mouseIDs = unique(T.mouse_id);
    for mouseIdx = 1:length(mouseIDs)
        currMouseID = mouseIDs(mouseIdx);
        currMouseAV = T(T.mouse_id == currMouseID,:);
        [~,id_idx,~] = intersect(currMouseAV.avIndex, avList);
        if ~isempty(id_idx)
            if ~exist(num2str(currMouseID),'dir')
                mkdir(num2str(currMouseID));
            end
            cd(num2str(currMouseID));
            trackIdx = unique(currMouseAV(id_idx,:).track_id);
            objectPoints = load(fullfile(pointListFolder, [num2str(currMouseID) '.mat']));
            objects = 1:size(objectPoints.pointList.pointList,1);
            % generate needed values
            bregma = allenCCFbregma(); % bregma position in reference data space
            atlas_resolution = 0.010; % mm
            matchedSlicesIdx = [];
            for object_num = 1:length(trackIdx)
                selected_object = trackIdx(object_num);
                % get the object points for the currently analyzed object    
                curr_objectPoints = objectPoints.pointList.pointList{selected_object,1}(:, [3 2 1]);
                curr_SliceIdxes = objectPoints.pointList.pointList{selected_object,2};
                % use the point's position in the atlas to get the AP, DV, and ML coordinates
                ap = -(curr_objectPoints(:,1)-bregma(1))*atlas_resolution;
                dv = (curr_objectPoints(:,2)-bregma(2))*atlas_resolution;
                ml = (curr_objectPoints(:,3)-bregma(3))*atlas_resolution;
                roi_location_curr = [ap dv ml];
                % initialize array of region annotations
                roi_annotation_curr = cell(size(curr_objectPoints,1),3); 
                % loop through every point to get ROI locations and region annotations
                for point = 1:size(curr_objectPoints,1)

                    % find the annotation, name, and acronym of the current ROI pixel
                    ann = av(curr_objectPoints(point,1),curr_objectPoints(point,2),curr_objectPoints(point,3));
                    name = st.safe_name{ann};
                    acr = st.acronym{ann};
                    if strcmp(acr,targetAeras{currAeraIdx})
                        matchedSlicesIdx = [matchedSlicesIdx;curr_SliceIdxes(point)];
                    end
                    roi_annotation_curr{point,1} = ann;
                    roi_annotation_curr{point,2} = name;
                    roi_annotation_curr{point,3} = acr;

                end
                roi_table = table(roi_annotation_curr(:,2),roi_annotation_curr(:,3), ...
                        roi_location_curr(:,1),roi_location_curr(:,2),roi_location_curr(:,3), roi_annotation_curr(:,1), ...
                    'VariableNames', {'name', 'acronym', 'AP_location', 'DV_location', 'ML_location', 'avIndex'});
                writetable(roi_table,['Probe ' num2str(selected_object) '.csv']);
            end
            matchedSlicesIdx = unique(matchedSlicesIdx);
            if isempty(matchedSlicesIdx)
                destinationFolder = pwd;
                targetResultFolder = [rtsFolderRoot '\' num2str(currMouseID)];
                copyfile(targetResultFolder,destinationFolder);
            else
                for sliceIdxx = 1:length(matchedSlicesIdx)
                    sliceNumber = matchedSlicesIdx(sliceIdxx);
                    % create Atlas viewer figure
                    f = figure('Name','Atlas Viewer'); 
                    f = atbnogui(f,currMouseID,av,st,tv,processedFolderRoot,sliceNumber);
                    saveas(f,fullfile(saveLocation,targetAeras{currAeraIdx},num2str(currMouseID),['slice ' num2str(sliceNumber) '.png']));
                    close all
                end
            end
            cd ..
        end
    end
    cd ..
end
        
                
                
                
        
        
    
    