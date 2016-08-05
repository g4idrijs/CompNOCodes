path('/home/tjh/git/zemp_lab/Matlab/bft',path);
path('/home/tjh/Matlab Simulator/',path);

load('/media/share/Tyler/Data/20130402-AnglesDopplerPA/thyroid1.mat');

f0 = Trans.frequency*1e6;              %  Central frequency                        [Hz]
fs = f0*Receive(1).samplesPerWave;     %  Sampling frequency                       [Hz]
c = Resource.Parameters.speedOfSound;              %  Speed of sound                           [m/s]
no_elements = Trans.numelements;      %  Number of elements in the transducer     

lambda = c / f0;       % Wavelength                                [m]
pitch = Trans.spacing*lambda;    % Pitch - center-to-center        [m]
width = Trans.elementWidth*lambda;     % Width of the element      [m]
kerf = pitch - width;  % Inter-element spacing                     [m]

%  Initialize the program

bft_init;

%  Set some paramters

bft_param('c', c);

bft_param('fs', fs);

xdc = bft_linear_array(no_elements, width, kerf);

%  Define and create the image
no_lines = 256;

imgWidth = PData(1).Size(2)*lambda;
xScanlines = linspace(-imgWidth/2, imgWidth/2, no_lines);

nSamples = 2*(Receive(1).endDepth-Receive(1).startDepth)*Receive(1).samplesPerWave;

if (exist('usTimeStart','var'))
    start_t = usTimeStart;
else
    start_t = 0;
end

bft_no_lines(no_lines);
for j = 1 : no_lines
    %bft_apodization(xdc, 0 , ones(1,no_elements),j);
    bft_center_focus([xScanlines(j) 0 0],j);
    bft_dynamic_focus(xdc, 0, 0,j);
end

if (~exist('na','var'))
    na=1;
end
fNum = 0;
for i=1:na
    rf_data = RcvData{1}(nSamples*(i-1)+1:nSamples*i,:,1);
    tic
    bf_temp = bft_beamform(start_t, rf_data);
    fprintf('bft: %f\n',toc);
    start_sample = floor(start_t*fs);
    if (start_sample==0)
        start_sample=1;
    end
    no_temp_samples = length(bf_temp);
    bf_lines(start_sample:start_sample+no_temp_samples-1,:) = bf_temp(1:no_temp_samples,:);
    env_bf = abs(hilbert(bf_lines));
    tic
    [x,y,img] = slidingDS(rf_data);
    fprintf('in-house: %f\n',toc);
    img = abs(hilbert(img));
    
    shifts = (xScanlines*sin(TX(i).Steer(1))/lambda)*Receive(1).samplesPerWave;
    shifts = shifts-min(shifts);
    for j=1:size(env_bf,2)
        env_bf(:,j) = interp1(1:size(env_bf,1),env_bf(:,j),(1:size(env_bf,1))+shifts(j));
        img(:,j) = interp1(1:size(img,1),img(:,j),(1:size(img,1))+shifts(j));
    end
    
    if (i==1)
        total_bf = env_bf;
    else
        total_bf = total_bf + env_bf;
    end
    
    if (i==1)
        total_img = img;
    else
        total_img = total_img + img;
    end
end

bft_end

%total_bf = total_bf/ max(max(abs(total_bf)));
figure;
imagesc(x*1000,y(1:1500)*1000,20*log10(total_bf(1:1500,:)+0.001));
set(gca,'Fontsize', 16);
ylabel('Depth (mm)','Fontsize', 16);
xlabel('Lateral position (mm)','Fontsize', 16);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'Visible','off');
title('Beamformed by BFT');
xLim = get(gca,'XLim'); yLim = get(gca,'YLim'); scale = (yLim(2)-yLim(1))/(xLim(2)-xLim(1));
Position = [0 0 3.25*2 3.25*scale*2];
set(gca,'Units','inches')
set(gca,'Position', Position);
set(gcf,'Units','inches');
set(gcf,'Position', Position);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition', Position);
set(gcf,'PaperSize', [Position(3) Position(4)]);
line([xLim(2)-6 xLim(2)-1], [yLim(1)+0.5 yLim(1)+0.5], 'LineWidth',4,'Color',[1 0 0]);
text(xLim(2)-3.5, yLim(1)+1, '5mm', 'FontSize', 32,'Color',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','top');
h=colorbar('West');
set(h,'YColor',[1 0 0]);
set(h,'XColor',[1 0 0]);
set(h,'FontSize',32);
colormap(gray);
print(gcf,'-dpng','bft_beamforming.png');

%total_img = total_img/ max(max(abs(total_img)));
figure;
imagesc(x*1000,y(1:1500)*1000,20*log10(total_img(1:1500,:)+0.001));
set(gca,'Fontsize', 16);
ylabel('Depth (mm)','Fontsize', 16);
xlabel('Lateral position (mm)','Fontsize', 16);
set(gca,'DataAspectRatio',[1 1 1]);
set(gca,'Visible','off');
title('Beamformed by Tyler');
xLim = get(gca,'XLim'); yLim = get(gca,'YLim'); scale = (yLim(2)-yLim(1))/(xLim(2)-xLim(1));
Position = [0 0 3.25*2 3.25*scale*2];
set(gca,'Units','inches')
set(gca,'Position', Position);
set(gcf,'Units','inches');
set(gcf,'Position', Position);
set(gcf,'PaperUnits','inches');
set(gcf,'PaperPosition', Position);
set(gcf,'PaperSize', [Position(3) Position(4)]);
line([xLim(2)-6 xLim(2)-1], [yLim(1)+0.5 yLim(1)+0.5], 'LineWidth',4,'Color',[1 0 0]);
text(xLim(2)-3.5, yLim(1)+1, '5mm', 'FontSize', 32,'Color',[1 0 0],'HorizontalAlignment','center','VerticalAlignment','top');
h=colorbar('West');
set(h,'YColor',[1 0 0]);
set(h,'XColor',[1 0 0]);
set(h,'FontSize',32);
colormap(gray);
print(gcf,'-dpng','tyler_beamforming_nofNum.png');