%drawTotalProbes
% directory of reference atlas files
annotation_volume_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\structure_tree_safe_2017.csv';
pointListFolder = 'H:\NP histology\probe_points\';

active_probe_length = 3.84;

% distance queried for confidence metric -- in um
probe_radius = 100; 

% overlay the distance between parent regions in gray (this takes a while)
show_parent_category = false; 

% plot this far or to the bottom of the brain, whichever is shorter -- in mm
distance_past_tip_to_plot = .5;

% set scaling e.g. based on lining up the ephys with the atlas
% set to *false* to get scaling automatically from the clicked points
scaling_factor = false;
% ---------------------
% additional parameters
% ---------------------
% plane used to view when points were clicked ('coronal' -- most common, 'sagittal', 'transverse')
plane = 'coronal';

% probe insertion direction 'down' (i.e. from the dorsal surface, downward -- most common!) 
% or 'up' (from a ventral surface, upward)
probe_insertion_direction = 'down';

% show a table of regions that the probe goes through, in the console
show_region_table = true;
      
% black brain?
black_brain = true;

fwireframe = [];

% scale active_probe_length appropriately
active_probe_length = active_probe_length*100;

ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0]; 

if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end

fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on; 
fwireframe.InvertHardcopy = 'off';

cd(pointListFolder);
pointLists = dir('*.mat');
if length(pointLists) > 11
    ProbeColors = .75*distinguishable_colors(length(pointLists),'k');
end

probe_lengths = cell(length(pointLists),2);
%% define the cell for probelengths here 
% eg£º
%probe_lengths{1,1} = 'probe_ponts_mouse01.mat'; % pointList for each mouse
%probe_lengths{1,2} = [1,2,3,4,5]; % corresponding implanting depth
% ...
%probe_lengths{length(pointLists),1} = ['probe_ponts_mouse' num2str(length(pointLists)) '.mat'];
%probe_lengths{length(pointLists),2} = [3,2.2,5,4.9,3.8,6];
%%

for i = 1:length(pointLists)
    currPointListName = load(pointLists(i).name);
    for j = 1:size(currPointListName.pointList.pointList,1)
        % get the probe points for the currently analyzed probe 
        if strcmp(plane,'coronal')
            currPointList = currPointListName.pointList.pointList{j,1}(:, [3 2 1]);
        elseif strcmp(plane,'sagittal')
            currPointList = currPointListName.pointList.pointList{j,1}(:, [1 2 3]);
        elseif strcmp(plane,'transverse')
            currPointList = currPointListName.pointList.pointList{j,1}(:, [1 3 2]);
        end
        currLengths = probe_lengths{string(cell2mat([probe_lengths(:,1)])) == pointLists(i).name,2};
        currLength = currLengths(j);
      % get the scaling-factor method to use
        if scaling_factor
            use_tip_to_get_reference_probe_length = false;
            reference_probe_length = currLength * scaling_factor;
            disp(['probe scaling of ' num2str(scaling_factor) ' determined by user input']);    
        else
            use_tip_to_get_reference_probe_length = true;
            disp(['getting probe scaling from histology data...']);
        end
        [m,p,s] = best_fit_line(currPointList(:,1), currPointList(:,2), currPointList(:,3));
        if isnan(m(1))
            disp(['no points found for probe ' num2str(j) ' in mouse ' num2str(i)]);
            continue
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
            m_further_up = m - p*50; % is there more brain 500 microns up along the track?
            ann_further_up = av(round(max(1,m_further_up(1))),round(max(1,m_further_up(2))),round(max(1,m_further_up(3))));
            if strcmp(st.safe_name(ann_further_up), 'root')
                out_of_brain = true;
            end
        end
    end

    % focus on wireframe plot
    figure(fwireframe);
    
    if use_tip_to_get_reference_probe_length
        % find length of probe in reference atlas space
        if strcmp(probe_insertion_direction, 'down')
            [depth, tip_index] = max(currPointList(:,2));
        elseif strcmp(probe_insertion_direction, 'up')
            [depth, tip_index] = min(currPointList(:,2));    
        end
        reference_probe_length_tip = sqrt(sum((currPointList(tip_index,:) - m).^2)); 
    
        % and the corresponding scaling factor
        shrinkage_factor = (reference_probe_length_tip / 100) / currLength;
    
        % display the scaling
        disp(['probe length of ' num2str(reference_probe_length_tip/100) ' mm in reference atlas space compared to a reported ' num2str(currLength) ' mm']);
        disp(['probe scaling of ' num2str(shrinkage_factor)]); disp(' ');
    
        % plot line the length of the probe in reference space
        probe_length_histo = round(reference_probe_length_tip);
    
        % if scaling_factor is user-defined as some number, use it to plot the length of the probe
    else 
        probe_length_histo = round(reference_probe_length * 100); 
    end

    % find the percent of the probe occupied by electrodes
    percent_of_tract_with_active_sites = min([active_probe_length / (currLength*100), 1.0]);
    active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);
    active_probe_position = round([active_site_start  probe_length_histo]);

    % plot line the length of the active probe sites in reference space
    plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
    'Color', ProbeColors(i,:), 'LineWidth', 1);   
    end
end