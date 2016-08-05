%BFT_ADD_IMAGE Add a low resolution to hi resolution image.
%
%USAGE : [hi_res] = bft_add_image(hi_res, lo_res, element, start_time)
%
%INPUT  : hi_res  - High resolution RF image. One column per scan line
%         lo_res  - Low resolution RF image. One column per scan line
%         element - Number of element, used to acquire the low resolution
%                   image
%         time    - Arrival time of the first sample of the RF lines.
%
%OUTPUT : hi_res - The high resolution image
%
%VERSION 1.0 Feb.21 2000, Svetoslav Nikolov

function hi_res1 = add_image(hi_res, lo_res, element, start_t)
hi_res1 = bft(13, hi_res, lo_res, element, start_t);
