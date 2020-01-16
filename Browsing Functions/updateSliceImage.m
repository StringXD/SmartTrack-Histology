function ud = updateSliceImage(ud,sliceNumber)

    title_ending = '';
    ud.deleteAll = 0;
    ud.slice_num = sliceNumber;
    processed_image_name = ud.processed_image_names{ud.slice_num};
    current_slice_image = flip(imread(fullfile(ud.processed_images_folder, processed_image_name)),1);
    if size(current_slice_image,1) > ud.ref_size(1)+2 || size(current_slice_image,2) > ud.ref_size(2)+2
        disp(['shrinking image to reference size ' num2str(ud.ref_size(1)) ' x ' num2str(ud.ref_size(2)) ' pxl'])
        current_slice_image = imresize(current_slice_image, ud.ref_size);
    end          
    set(ud.im, 'CData', current_slice_image); 


    file_transformations = fullfile(ud.processed_images_folder, 'transformations\\' ,...
                            [processed_image_name(1:end-4) '_transform_data.mat']);

    set(ud.pointHands(:), 'Visible', 'off'); 

    if exist(file_transformations,'file')
        % load transform data
        transform_data = load(file_transformations);
        transform_data = transform_data.save_transform;
        if ~isempty(transform_data.transform_points{2})
            ud.pointList = transform_data.transform_points{2};
            title_ending = ' (transform points loaded)';
        end
    else
        ud.pointList = [];
    end
    title(['Slice Viewer -- Slice ' num2str(ud.slice_num) '/' num2str(ud.total_num_files) title_ending])    