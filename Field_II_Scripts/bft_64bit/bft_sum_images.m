%BFT_SUM_IMAGES Sum 2 low resolution images in 1 high resolution.
%
%USAGE  : [hi_res] = bft_sum_images(image1, ele1, image2, ele2, time)
%
%INPUT  : image1 - Matrix with the RF data for the image. The number
%                  of columns corresponds to the number of lines
%         ele1   - Number of emitting element used to obtain the image. 
%         image2 - Matrix with the RF data for the image. The number
%                  of columns corresponding to the number of lines
%         ele2   - Number of emitting element used to obtain the image.
%         time   -  The arrival time of the first samples. The two images
%                   must be aligned in time
%
%OUTPUT : hi_res - Higher resolution image
%
%VERSION:1.0 Feb.18 2000, Svetoslav Nikolov
function [hi_res] = bft_sum_images(image1, ele1, image2, ele2, time)
[hi_res] = bft(12, image1, ele1, image2, ele2, time);
