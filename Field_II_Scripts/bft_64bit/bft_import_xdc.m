%BFT_IMPORT_XDC Import transducer definition from Field II
%
%USAGE   : xdc_bft = bft_import_xdc(xdc_field)
%
%INPUTS  : xdc_field - Transducer definition from Field II
% 
%OUTPUT  : xdc_bft - Transducer definition in BFT
%
%CREATED : 06 Oct 2000, Svetoslav Nikolov
%

function xdc_bft = bft_import_xdc(xdc_field)

data = xdc_get(xdc_field, 'rect');

x = data(24,:);
y = data(25,:);
z = data(26,:);

xdc_bft = bft_transducer([x' y' z']);

