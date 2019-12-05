function[] = Display_Designed_Track(matrix)
%% the structure of matrix:
% Each row: AP ML DV angle
%% ENTER PARAMETERS AND FILE LOCATION

% directory of reference atlas files
annotation_volume_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\annotation_volume_10um_by_index.npy';
structure_tree_location = 'E:\prJ\neuropixels\histology location analysis\allenCCF\structure_tree_safe_2017.csv';

% directory of saving result figures and tables
save_location = 'E:\test\fig';

% --------------
% key parameters
% --------------
% how far into the brain did you go from the surface, either for each probe or just one number for all -- in mm
probe_lengths = matrix(:,3)';

% from the bottom tip, how much of the probe contained recording sites -- in mm
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

% probe insertion direction 'down' (i.e. from the dorsal surface, downward -- most common!) 
% or 'up' (from a ventral surface, upward)
probe_insertion_direction = 'down';

% show a table of regions that the probe goes through, in the console
show_region_table = true;
      
% black brain?
black_brain = true;

%% GET AND PLOT PROBE VECTOR IN ATLAS SPACE

% load the reference brain annotations
if ~exist('av','var') || ~exist('st','var')
    disp('loading reference atlas...')
    av = readNPY(annotation_volume_location);
    st = loadStructureTree(structure_tree_location);
end
% Color design
if size(matrix,1) > 11
    ProbeColors = .75*distinguishable_colors(size(matrix,1),'k');
else
    ProbeColors = .75*[1.3 1.3 1.3; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 .9; .65 .4 .25; .7 .95 .3; .7 0 0; .6 0 .7; 1 .6 0];
end

% order of colors: {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','purple','orange','red'};
fwireframe = [];

% scale active_probe_length appropriately
active_probe_length = active_probe_length*100;

%% PLOT EACH PROBE -- FIRST FIND ITS TRAJECTORY IN REFERENCE SPACE

% create a new figure with wireframe
fwireframe = plotBrainGrid([], [], fwireframe, black_brain);
hold on; 
fwireframe.InvertHardcopy = 'off';

% define bregma and scale
bregma = [540,570,0]; % AP LR DV
atlas_res = 0.010;
%atlas_res = [0.01,4.25/475,0.01];

for selected_probe = 1:size(matrix,1)
    
    % get user-defined probe length from experiment
    if length(probe_lengths) > 1
        probe_length = probe_lengths(selected_probe);
    else
        probe_length = probe_lengths;
    end
    
    % determine "origin" at top of brain -- step upwards along tract direction until tip of brain / past cortex
    % m is 2 mm below bregma height; p is the eigenvector for largest eigenvalue
    matrix(selected_probe,1) = -matrix(selected_probe,1);
    m = bregma + matrix(selected_probe,1:3)/atlas_res;
    m = [m(1),m(3),m(2)]; % AP DV LR
    p = [0,1,0];
    
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
    
    % compute p vector based on implantation angle
    p = [0,cosd(matrix(selected_probe,4)),-sind(matrix(selected_probe,4))];
    
    % focus on wireframe plot
    figure(fwireframe);
    
    % plot brain entry point
    plot3(m(1), m(3), m(2), 'r*','linewidth',1)
    plot3(m(1), 1140-m(3), m(2), 'r*','linewidth',1)
    
    probe_length_histo = round(probe_length * 100);
    
    % find the percent of the probe occupied by electrodes
    percent_of_tract_with_active_sites = min([active_probe_length / (probe_length*100), 1.0]);
    active_site_start = probe_length_histo*(1-percent_of_tract_with_active_sites);
    active_probe_position = round([active_site_start  probe_length_histo]);
    
    % plot line the length of the active probe sites in reference space
    plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], m(3)+p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);
    plot3(m(1)+p(1)*[active_probe_position(1) active_probe_position(2)], 1140-m(3)-p(3)*[active_probe_position(1) active_probe_position(2)], m(2)+p(2)*[active_probe_position(1) active_probe_position(2)], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1);
    % plot line the length of the entire probe in reference space
    plot3(m(1)+p(1)*[1 probe_length_histo], m(3)+p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');
    plot3(m(1)+p(1)*[1 probe_length_histo], 1140-m(3)-p(3)*[1 probe_length_histo], m(2)+p(2)*[1 probe_length_histo], ...
        'Color', ProbeColors(selected_probe,:), 'LineWidth', 1, 'LineStyle',':');
    
    %% Get and plot brain region labels along the extent of each probe
    
    % convert error radius into mm
    error_length = round(probe_radius / 10);
    % find and regions the probe goes through, confidence in those regions, and plot them
    borders_table = plotDistToNearestToTip(m, p, av, st, probe_length_histo, error_length, active_site_start, distance_past_tip_to_plot, show_parent_category, show_region_table); % plots confidence score based on distance to nearest region along probe
    writetable(borders_table,[save_location '\Probe ' num2str(selected_probe) '.csv']);
    title(['Probe ' num2str(selected_probe)],'color',ProbeColors(selected_probe,:))
    saveas(gcf,[save_location '\Probe ' num2str(selected_probe) '.png']);
    close(gcf);


    pause(.05)
end
saveas(fwireframe,[save_location '\overview.fig']);
