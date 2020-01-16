function updateSliceNoGUI(f,slice_num,allData, slice_figure, save_location)

ud = get(f, 'UserData');

% scroll through slices
if ud.scrollMode == 3
    set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
    ud.showOverlay = 0;
    delete(ud.overlayAx); ud.overlayAx = [];  
    ud_slice = get(slice_figure, 'UserData');
    ud_slice.pointList = [];
    ud.slice_shift = ud.slice_shift + slice_num;
    slice_name = ud_slice.processed_image_names{ud.slice_at_shift_start+ud.slice_shift}(1:end-4);
    folder_transformations = fullfile(save_location, ['transformations' filesep]);
    for probe = 1:size(ud.pointList,1)
        set(ud.pointHands{probe, 3}(:),'Visible','off'); ud.pointHands{probe, 3} = [];
        for probe_point = 1:size(ud.pointList{probe,1},1)
            slice_point_belongs_to = ud.pointList{probe, 2}(probe_point);
            if slice_point_belongs_to == ud.slice_at_shift_start+ud.slice_shift
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'on');
            else
                set(ud.pointHands{probe, 1}(probe_point), 'Visible', 'off');
            end
        end
    end 
    ud.current_pointList_for_transform = zeros(0,2);
    try
        load([folder_transformations slice_name '_transform_data.mat']);
        % load transform data
        transform_data = load(fullfile(folder_transformations, [slice_name '_transform_data.mat']));  
        transform_data = transform_data.save_transform;
        % load new transform
        ud.transform = transform_data.transform;
        
        if ~isempty(transform_data.transform_points{1}) && ~isempty(transform_data.transform_points{2})
            ud.current_pointList_for_transform = transform_data.transform_points{1};
            ud_slice.pointList = transform_data.transform_points{2};
        else
            ud_slice.pointList = [];           
        end
        set(slice_figure, 'UserData', ud_slice);
        % load allen ref location
        ud.currentSlice = transform_data.allen_location{1}; ud.currentAngle = transform_data.allen_location{2};

        % create transformed histology image
        ud.curr_slice_trans = imread([folder_transformations slice_name '_transformed.tif']);
       
        % update figure
        set(f, 'UserData', ud);
        ud.loaded = true;
        
        ud.curr_slice_num = ud.slice_at_shift_start+ud.slice_shift;
        
        ud.histology_overlay = 1;
        
        set(ud.text,'Visible','on');
        fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)],'color','white');
    catch
        % if no transform, just show reference
        ud.histology_overlay = 0;
        ud.current_pointList_for_transform = zeros(0,2);
        set(ud.im, 'CData', ud.ref);
        ud.curr_im = ud.ref; set(f, 'UserData', ud);   
        set(ud.text,'Visible','off');
        fill([5 5 250 250],[5 50 50 5],[0 0 0]); ud.text(end+1) = text(5,15,['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift) ' - no transform'],'color','white');
    end
end
% update coordinates at the top
pixel = getPixel(ud.atlasAx);
updateStereotaxCoords(f,ud.currentSlice, pixel, ud.bregma, ud.bregmaText, ud.angleText, ud.currentSlice, ud.currentAngle(1), ud.currentAngle(2), ud.ref_size);
% if no angle, just change reference slice
if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
    reference_slice = squeeze(allData.tv(ud.currentSlice,:,:));
    ud.im_annotation = squeeze(allData.av(ud.currentSlice,:,:));   
    if ud.viewColorAtlas
        set(ud.im, 'CData', ud.im_annotation);
    else
        set(ud.im, 'CData', reference_slice);
    end
   % update title/overlay with brain region
    [name, acr, ann] = getPixelAnnotation(allData, pixel, ud.currentSlice);
    updateTitle(ud.atlasAx, name, acr)
    if ud.showOverlay    
        updateOverlay(f, allData, ann, slice_figure, save_location);
    end  
    ud.ref = uint8(reference_slice);
    set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
    ud.offset_map = zeros(ud.ref_size); 
    set(f, 'UserData', ud);
else
    image_size = size(squeeze(allData.av(ud.currentSlice,:,:)));
    angle_slice = zeros(image_size);
    if ud.currentAngle(1)==0; offset_DV = 0;
    else; offset_DV = -ud.currentAngle(1):sign(ud.currentAngle(1)):ud.currentAngle(1);
    end; start_index_DV = 1; 
    % loop through AP offsets
    num_DV_iters_add_ind = floor( (image_size(1) - floor( image_size(1) / length(offset_DV))*length(offset_DV)) / 2);
    for curr_DV_iter = 1:length(offset_DV)
      cur_offset_DV = offset_DV(curr_DV_iter);
      if cur_offset_DV == ud.currentAngle(1)
          end_index_DV = image_size(1);
      elseif curr_DV_iter <= num_DV_iters_add_ind  || length(offset_DV - curr_DV_iter) < num_DV_iters_add_ind
          end_index_DV = start_index_DV + floor( image_size(1) / length(offset_DV));
      else
          end_index_DV = start_index_DV + floor( image_size(1) / length(offset_DV)) - 1;
      end
      
       if ud.currentAngle(2)==0;  offset_ML = 0;
       else; offset_ML = -ud.currentAngle(2):sign(ud.currentAngle(2)):ud.currentAngle(2);
       end; start_index_ML = 1;
    % nested: loop through ML offsets
    num_ML_iters_add_ind = floor( (image_size(2) - floor( image_size(2) / length(offset_ML))*length(offset_ML)) / 2);
    for curr_ML_iter = 1:length(offset_ML)
        cur_offset_ML = offset_ML(curr_ML_iter);
        if cur_offset_ML == ud.currentAngle(2)
            end_index_ML = image_size(2);
        elseif curr_ML_iter <= num_ML_iters_add_ind  || length(offset_ML - curr_ML_iter) < num_ML_iters_add_ind
            end_index_ML = start_index_ML + floor( image_size(2) / length(offset_ML));
        else
            end_index_ML = start_index_ML + floor( image_size(2) / length(offset_ML)) - 1;
        end
        % update current slice
        try
        angle_slice(start_index_DV:end_index_DV, start_index_ML:end_index_ML) = ...
         squeeze(allData.tv(ud.currentSlice + cur_offset_DV + cur_offset_ML,start_index_DV:end_index_DV,start_index_ML:end_index_ML));
        catch
            disp('')
        end
        ud.im_annotation(start_index_DV:end_index_DV,start_index_ML:end_index_ML) = squeeze(allData.av(ud.currentSlice + cur_offset_DV + cur_offset_ML,...
                                                            start_index_DV:end_index_DV,start_index_ML:end_index_ML));
        ud.offset_map(start_index_DV:end_index_DV, start_index_ML:end_index_ML) = cur_offset_DV + cur_offset_ML;
        start_index_ML = end_index_ML + 1;
    end
        start_index_DV = end_index_DV + 1;
    end
    if ud.viewColorAtlas
        set(ud.im, 'CData', ud.im_annotation);
    elseif ~ud.showAtlas  
        set(ud.im, 'CData', angle_slice);
    end
    ud.ref = uint8(angle_slice);
  set(ud.pointHands_for_transform(:), 'Visible', 'off'); 
end
% in all cases. update histology overlay
if ud.histology_overlay == 1 || ud.histology_overlay == 2
    updateHistology(f,ud); ud = get(f, 'UserData');
else
    if ud.viewColorAtlas
        ud.curr_im = ud.im_annotation;
    else
        ud.curr_im = ud.ref;
    end    
end

% then update boundary overlay
if ud.showAtlas
    updateBoundaries(f,ud, allData)
end
set(f, 'UserData', ud);


 
  
function updateHistology(f, ud)
    if ud.histology_overlay == 2
        image_blend =  imfuse(uint8(ud.ref*.6), ud.curr_slice_trans(:,:,:),'blend','Scaling','none');
        set(ud.im, 'CData', image_blend);
        ud.curr_im = image_blend;
    elseif ud.histology_overlay == 1
        set(ud.im, 'CData', ud.curr_slice_trans);
        ud.curr_im = ud.curr_slice_trans;
    end
    set(f, 'UserData', ud);
    
function updateBoundaries(f, ud, allData)
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


function updateStereotaxCoords(f,currentSlice, pixel, bregma, bregmaText, angleText, slice_num, ap_angle, ml_angle, ref_size)
ud = get(f, 'UserData');
set(bregmaText,'String',['Slice ' num2str(ud.slice_at_shift_start+ud.slice_shift)]);
set(angleText, 'String', ['Slice ' num2str((bregma(1) - slice_num)/100) ' AP, DV angle ' num2str(round(atand(ap_angle/(ref_size(1)/2)),1)) '^{\circ}, ML angle ' num2str(round(atand(ml_angle/570),1)) '^{\circ}']);

function pixel = getPixel(ax)

currPoint = get(ax,'currentpoint');  % The current point w.r.t the axis.

Cx = currPoint(1,1); Cy = currPoint(1,2);
pixel = round([Cy Cx]);

function updateOverlay(f, allData, ann, slice_figure, save_location)
ud = get(f, 'UserData');
if isempty(ud.overlayAx) % first time
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        avo = plotAVoverlay(fliplr(squeeze(allData.av(ud.currentSlice,:,:))'), ann, ud.atlasAx);
    else
        avo = plotAVoverlay(fliplr(ud.im_annotation)', ann, ud.atlasAx);
    end
    ud.overlayAx = avo;
    set(ud.overlayAx, 'HitTest', 'off');
    set(f, 'UserData', ud);
else
    ovIm = get(ud.overlayAx, 'Children');
%     set(ovIm, 'HitTest', 'off');
    
    if ud.currentAngle(1) == 0 && ud.currentAngle(2) == 0
        thisSlice = squeeze(allData.av(ud.currentSlice,:,:));
    else
        thisSlice = ud.im_annotation;
    end

    set(ovIm, 'CData', flipud(thisSlice));    
    plotAVoverlay(fliplr(thisSlice'), ann, ud.atlasAx, ud.overlayAx);
    set(ovIm, 'ButtonDownFcn', @(f,k)atlasClickCallback(f, k, slice_figure, save_location));

end

function [name, acr, ann] = getPixelAnnotation(allData, pixel, currentSlice)
if pixel(1)>0&&pixel(1)<size(allData.av,2) && pixel(2)>0&&pixel(2)<=size(allData.av,3)
    ann = allData.av(currentSlice,pixel(1),pixel(2));
    name = allData.st.safe_name(ann);
    acr = allData.st.acronym(ann);
else
    ann = []; name = []; acr = [];
end

function updateTitle(ax, name, acr)
if ~isempty(name)
    title(ax, [name{1} ' (' acr{1} ')']);
else
    title(ax, 'not found');
end
