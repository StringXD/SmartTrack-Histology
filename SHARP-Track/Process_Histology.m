% ------------------------------------------------------------------------
%          Crop, Rotate, Adjust Contrast, and Downsample Histology
% ------------------------------------------------------------------------


%%  SET FILE AND PARAMETERS

% * remember to run one cell at a time, instead of the whole script at once *

% directory of histology images
Image4register_folder = 'E:\prJ\neuropixels\NP histology\xd_20191120\25\reg';
Image4ROIs_folder = 'E:\prJ\neuropixels\NP histology\xd_20191120\25\roi';
% directory to save the processed images -- can be the same as the above image_folder
% results will be put inside a new folder called 'processed' inside of this image_folder
save_folder = 'E:\prJ\neuropixels\NP histology\xd_20191120\25';

% name of images, in order anterior to posterior or vice versa
% once these are downsampled they will be named ['original name' '_processed.tif']
image4register_file_names = dir([Image4register_folder filesep '*.tif']); % get the contents of the image_folder
image4register_file_names = natsortfiles({image4register_file_names.name});
image4rois_file_names = dir([Image4ROIs_folder filesep '*.tif']); % get the contents of the image_folder
image4rois_file_names = natsortfiles({image4rois_file_names.name});
% image_file_names = {'slide no 2_RGB.tif','slide no 3_RGB.tif','slide no 4_RGB.tif'}; % alternatively, list each image in order

% if the images are individual slices (as opposed to images of multiple
% slices, which must be cropped using the cell CROP AND SAVE SLICES)
image_files_are_individual_slices = true;

% use images that are already at reference atlas resolution (here, 10um/pixel)
use_already_downsampled_image = false; 

% pixel size parameters: microns_per_pixel of large images in the image
% folder (if use_already_downsampled_images is set to false);
% microns_per_pixel_after_downsampling should typically be set to 10 to match the atlas
microns_per_pixel = 3.25;
microns_per_pixel_after_downsampling = 10;


% ----------------------
% additional parameters
% ----------------------

% if the images are cropped (image_file_are_individual_slices = false),
% name to save cropped slices as; e.g. the third cropped slice from the 2nd
% image containing many slices will be saved as: save_folder/processed/save_file_name02_003.tif
save_file_name = 'M25';

% increase gain if for some reason the images are not bright enough
gain = 1; 

% size in pixels of reference atlas brain coronal slice, typically 800 x 1140
atlas_reference_size = [800 1140]; 







% finds or creates a folder location for processed images -- 
% a folder within save_folder called processed
folder_processed_images4reg = fullfile(save_folder, 'processed_4reg');
if ~exist(folder_processed_images4reg)
    mkdir(folder_processed_images4reg)
end
folder_processed_images4roi = fullfile(save_folder, 'processed_roi');
if ~exist(folder_processed_images4roi)
    mkdir(folder_processed_images4roi)
end

%% LOAD AND PROCESS SLICE PLATE IMAGES

% close all figures
close all
   

% if the images need to be downsampled to 10um pixels (use_already_downsampled_image = false), 
% this will downsample and allow you to adjust contrast of each channel of each image from image_file_names
%
% if the images are already downsampled (use_already_downsampled_image = true), this will allow
% you to adjust the contrast of each channel
%
% Open Histology Viewer figure
try; figure(histology_figure);
catch; histology_figure = figure('Name','Histology Viewer'); end
warning('off', 'images:initSize:adjustingMag'); warning('off', 'MATLAB:colon:nonIntegerIndex');

% Function to downsample and adjust histology image
HistologyBrowser(histology_figure, save_folder, Image4register_folder, image4register_file_names, folder_processed_images4reg, image_files_are_individual_slices, ...
use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)
if ~strcmp(Image4register_folder,Image4ROIs_folder)    
    pause
    close all
    histology_figure = figure('Name','Histology Viewer');
    HistologyBrowser(histology_figure, save_folder, Image4ROIs_folder, image4rois_file_names, folder_processed_images4roi, image_files_are_individual_slices, ...
            use_already_downsampled_image, microns_per_pixel, microns_per_pixel_after_downsampling, gain)
end
  

%% CROP AND SAVE SLICES -- run once the above is done, if image_file_are_individual_slices = false

% close all figures
close all

image_file_names = image4register_file_names;


% run this function if the images from image_file_names have several
% slices, per image (e.g. an image of an entire histology slide)
%
% this function allows you to crop all the slices you would like to 
% process im each image, by drawing rectangles around them in the figure. 
% these will be further processed in the next cell
if ~image_files_are_individual_slices
    histology_figure = figure('Name','Histology Viewer');
    HistologyCropper(histology_figure, save_folder, image_file_names, atlas_reference_size, save_file_name, use_already_downsampled_image)
else
    disp('individually cropped slices already available')
end


%% GO THROUGH TO FLIP HORIZONTAL SLICE ORIENTATION, ROTATE, SHARPEN, and CHANGE ORDER

% close all figures
close all
            
% this takes images from folder_processed_images ([save_folder/processed]),
% and allows you to rotate, flip, sharpen, crop, and switch their order, so they
% are in anterior->posterior or posterior->anterior order, and aesthetically pleasing
% 
% it also pads images smaller than the reference_size and requests that you
% crop images larger than this size
%
% note -- presssing left or right arrow saves the modified image, so be
% sure to do this even after modifying the last slice in the folder
if ~strcmp(Image4register_folder,Image4ROIs_folder)
    register_figure = figure('Name','Slice Viewer');
    SliceFlipper(register_figure, folder_processed_images4roi, atlas_reference_size,1,folder_processed_images4reg)
else
    register_figure = figure('Name','Slice Viewer');
    SliceFlipper(register_figure, folder_processed_images4reg, atlas_reference_size,0)
end


