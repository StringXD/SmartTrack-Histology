function f = atbnogui(f,mouse_id,av,st,tv,processedFolderRoot,sliceNumber)
% directory of histology
processed_images_folder = fullfile(processedFolderRoot,num2str(mouse_id),'processed',filesep);  
processed_roi_folder = fullfile(processedFolderRoot,num2str(mouse_id),'processed',filesep);
% name the saved probe points, to avoid overwriting another set of probes going in the same folder
probe_save_name_suffix = 'electrode_1'; 

% plane to view ('coronal', 'sagittal', 'transverse')
plane = 'coronal';

% select the plane for the viewer
if strcmp(plane,'coronal')
    av_plot = av;
    tv_plot = tv;
elseif strcmp(plane,'sagittal')
    av_plot = permute(av,[3 2 1]);
    tv_plot = permute(tv,[3 2 1]);
elseif strcmp(plane,'transverse')
    av_plot = permute(av,[2 3 1]);
    tv_plot = permute(tv,[2 3 1]);
end

% show histology in Slice Viewer
try; figure(slice_figure_browser); title('');
catch; slice_figure_browser = figure('Name','Slice Viewer'); end
reference_size = size(tv);
% slice browser begins here
slice_figure = slice_figure_browser;
% initialize user data variables held by the figure
processed_images = dir([processed_images_folder filesep '*.tif']);
ud_slice.processed_image_names = natsortfiles({processed_images.name});
total_num_files = size(processed_images,1); disp(['found ' num2str(total_num_files) ' processed slice images']);
ud_slice.total_num_files = total_num_files;
ud_slice.break = 0; 
ud_slice.slice_num = 1;
ud_slice.key = 1; 
ud_slice.pointList = []; 
ud_slice.pointHands = [];
ud_slice.getPoint = 0;
ud_slice.ref_size = reference_size(2:3);

ud_slice.processed_images = processed_images;
ud_slice.processed_images_folder = processed_images_folder;
ud_slice.sliceAx = axes('Position', [0.05 0.05 0.9 0.9]);
hold(ud_slice.sliceAx, 'on');
set(ud_slice.sliceAx, 'HitTest', 'off');
ud_slice.im = plotTVslice(zeros(ud_slice.ref_size, 'uint8'));

try; screen_size = get(0, 'ScreenSize'); screen_size = screen_size(1,3:4)./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end
set(slice_figure,'Position', [150*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
movegui(slice_figure,'onscreen')
% atlas browser starts here
templateVolume = tv_plot;
annotationVolume = av_plot;
structureTree = st;
slice_figure = slice_figure_browser;
save_location = processed_roi_folder;
save_suffix = probe_save_name_suffix;
%
figure(f);
try; screen_size = get(0, 'ScreenSize'); screen_size = screen_size(1,3:4)./[2560 1440];
catch; screen_size = [1900 1080]./[2560 1440];
end 
set(f,'Position', [1050*screen_size(1) 660*screen_size(2) 880*screen_size(1) 650*screen_size(2)])
movegui(f,'onscreen')
% initialize user data variables held by the figure
ud.bregma = allenCCFbregma; 
ud.currentSlice = ud.bregma(1); 
ud.currentAngle = zeros(2,1);
ud.scrollMode = 3;
ud.transform_type = 'projective'; %can change to 'affine' or 'pwl'
ud.oldContour = [];
ud.showContour = false;
ud.showOverlay = false; ud.overlayAx = [];
ud.pointList = cell(1,3); ud.pointList{1} = zeros(0,3); 
ud.pointHands = cell(1,3);
ud.probe_view_mode = false;
ud.currentProbe = 0; ud.ProbeColors = [1 1 1; 1 .75 0;  .3 1 1; .4 .6 .2; 1 .35 .65; .7 .7 1; .65 .4 .25; .7 .95 .3; .7 0 0; .5 0 .6; 1 .6 0]; 
ud.ProbeColors = [ud.ProbeColors; distinguishable_colors(5,ud.ProbeColors)];
ud.ProbeColor =  {'white','gold','turquoise','fern','bubble gum','overcast sky','rawhide', 'green apple','red','purple','orange', 'ex1','ex2','ex3','ex4','ex5'};
ud.getPoint_for_transform =false; ud.pointList_for_transform = zeros(0,2); ud.pointHands_for_transform = [];
ud.current_pointList_for_transform = zeros(0,2); ud.curr_slice_num = 1;
ud.clicked = false;
ud.showAtlas = true;
ud.viewColorAtlas = false;
ud.histology_overlay = 0; 
ud.atlasAx = axes('Position', [0.05 0.05 0.9 0.9]);
ud.transform = [];
ud.transformed_slice_figure = [];
ud.slice_shift = 0;
ud.loaded_slice = 0;
ud.slice_at_shift_start = 1;
ud.text = [];
ud.scrollBackwards = false;
reference_image = squeeze(templateVolume(ud.currentSlice,:,:));
ud.im = plotTVslice(reference_image);
ud.ref_size = size(reference_image);
ud.ref = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_im = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.curr_slice_trans = uint8(squeeze(templateVolume(ud.currentSlice,:,:)));
ud.im_annotation = squeeze(annotationVolume(ud.currentSlice,:,:));
ud.atlas_boundaries = zeros(ud.ref_size,'uint16');
ud.offset_map = zeros(ud.ref_size);
ud.loaded = 0;
ud.bregmaText = annotation('textbox', [0 0.95 0.4 0.05], ...
    'String', '[coords]', 'EdgeColor', 'none', 'Color', 'k');

ud.angleText = annotation('textbox', [.7 0.95 0.4 0.05], ...
    'EdgeColor', 'none', 'Color', 'k');

allData.tv = templateVolume;
allData.av = annotationVolume;
allData.st = structureTree;
hold(ud.atlasAx, 'on');
set(ud.atlasAx, 'HitTest', 'off');
set(f, 'UserData', ud);
% retrieve user data from figure
ud = get(f, 'UserData');
% update slice_num in slice browser
ud_slice = updateSliceImage(ud_slice,sliceNumber); 
set(slice_figure,'UserData',ud_slice);
% extract data of certain slice
ud_slice = get(slice_figure, 'UserData');
updateSliceNoGUI(f, sliceNumber - 1, allData, slice_figure, save_location);ud = get(f, 'UserData');
slice_name = ud_slice.processed_image_names{ud_slice.slice_num}(1:end-4);
folder_transformations = fullfile(save_location, ['transformations' filesep]);
try
        ud.curr_slice_num = ud_slice.slice_num;
        % remove overlay
        ud.showOverlay = 0;
        delete(ud.overlayAx); ud.overlayAx = [];
        % load transform data
        transform_data = load(fullfile(folder_transformations, [slice_name '_transform_data.mat'])); 
        transform_data = transform_data.save_transform;
        % load new transform
        ud.transform = transform_data.transform;
        if ~isempty(transform_data.transform_points{1}) && ~isempty(transform_data.transform_points{2})
            ud.current_pointList_for_transform = transform_data.transform_points{1};
            ud_slice.pointList = transform_data.transform_points{2};
            set(slice_figure, 'UserData', ud_slice);
        end
        % load allen ref location
        ud.currentSlice = transform_data.allen_location{1}; ud.currentAngle = transform_data.allen_location{2};
        % create transformed histology image
        current_slice_image = flip(get(ud_slice.im, 'CData'),1); R = imref2d(size(ud.ref));
        ud.curr_slice_trans = imwarp(current_slice_image, ud.transform, 'OutputView',R);
        % update figure
        temp_scroll_mode = ud.scrollMode; ud.scrollMode = 4; set(f, 'UserData', ud);
        updateSliceNoGUI(f, 0, allData, slice_figure, save_location);ud = get(f, 'UserData');
        ud.scrollMode = temp_scroll_mode;
        ud.loaded = true;
        ud.slice_at_shift_start = ud_slice.slice_num; ud.slice_shift = 0;
        ud.loaded_slice = ud_slice.slice_num;            
     % load probe points
        if ~size(ud.pointList{1,1},1)
            probe_points = load(fullfile(save_location, ['probe_points' save_suffix]));  disp('probe points loaded')
            ud.pointList = probe_points.pointList.pointList;
            ud.pointHands = probe_points.pointList.pointHands;
        else
            disp('probe points not loaded -- there are already some current probe points')
        end
    for probe = 1:size(ud.pointList,1)
        % create point plot handles anew
        try; set(ud.pointHands{probe, 1}(:),'Visible','off'); end
        for probe_point = 1:size(ud.pointList{probe, 1},1)
            ud.pointHands{probe, 1}(probe_point) = scatter(ud.atlasAx, ...
                ud.pointList{probe,1}(probe_point,1), ud.pointList{probe,1}(probe_point,2), 20, 'ro', ...
            'MarkerFaceColor', [.1 .1 .1],'MarkerEdgeColor', ud.ProbeColors(probe, :), ...
            'LineWidth',2);
        % set the point plot from the current slice visible
        slice_point_belongs_to = ud.pointHands{probe, 2}(probe_point);
        if slice_point_belongs_to == ud_slice.slice_num
            set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'on');
        else
            set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'off'); 
        end
        end
        % turn off best fit line if applicable
        set(ud.pointHands{probe, 3}(:),'Visible','off'); ud.pointHands{probe, 3} = [];
    end
    ud.slice_shift = 0;
catch
    disp(['loading failed']);
end
% add boundaries
if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
    curr_annotation = squeeze(allData.av(ud.currentSlice,:,:));
else
    curr_annotation = ud.im_annotation;
end
    
    atlas_vert_1 = double(curr_annotation(1:end-2,:));
    atlas_vert_2 = double(curr_annotation(3:end,:));
    atlas_vert_offset = abs( atlas_vert_1 - atlas_vert_2 ) > 0;
    shifted_atlas_vert1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert1(3:end,:) = atlas_vert_offset;
    shifted_atlas_vert2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_vert2(1:end-2,:) = atlas_vert_offset;

    atlas_horz_1 = double(curr_annotation(:,1:end-2));
    atlas_horz_2 = double(curr_annotation(:,3:end));
    atlas_horz_offset = abs( atlas_horz_1 - atlas_horz_2 )>0;
    shifted_atlas_horz1 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz1(:,3:end) = atlas_horz_offset;
    shifted_atlas_horz2 = zeros(size(curr_annotation(:,:)));
    shifted_atlas_horz2(:,1:end-2) = atlas_horz_offset;

    shifted_atlas = shifted_atlas_horz1 + shifted_atlas_horz2 + shifted_atlas_vert1 + shifted_atlas_vert2;

    atlas_boundaries = (shifted_atlas>0); ud.atlas_boundaries = atlas_boundaries;

    if ud.showAtlas
        image_blend =  uint8( imfuse(ud.curr_im, atlas_boundaries/3.5*(1+.35*isa(ud.curr_im,'uint16')),'blend','Scaling','none') )* 2;
        set(ud.im, 'CData', image_blend); 
    end
    
    set(f, 'UserData', ud);



