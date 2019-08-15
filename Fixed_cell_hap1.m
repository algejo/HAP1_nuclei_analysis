% Image extraction
function Fixed_cell_hap1(row,col,site,day)

% row = 4;col = 4;site = 1;day = 1;
tic
imagepath = 'C:\Users\Michael\';
shadingpath='C:\Users\Michael\Documents\MATLAB\Scripts\Fixed Cell Scripts\ShadingImages\DCYTC 10x\';
experimentpath='Documents\MATLAB\Experiments\ddx21\'; 
datadir = [imagepath,experimentpath,'Data\'];
rawdir = [imagepath, experimentpath,'Raw\'];
biasdir=[rawdir,'Bias\'];

% imagepath = '/mnt/fluffy/michael/';
% shadingpath='/mnt/fluffy/michael/Scripts/ShadingImages/DCYTC 10x/';
% experimentpath='/Fixed Cell Experiments/2017/20170325-shRNA-p50/'; 
% datadir = [imagepath,experimentpath,'Data/'];
% rawdir = [imagepath, experimentpath,'Raw/'];
% biasdir=[rawdir,'Bias/'];

% 'D:\ShadingImages\DCYTC 10x\CameraNoise_1080_bin2.mat';
% load([shadingpath,'CameraNoise_2160_bin1.mat'],'BG'); bgcmos = BG;
load([shadingpath,'CameraNoise_1080_bin2.mat'],'BG'); bgcmos = BG;
if ~exist(datadir,'dir')
    mkdir(datadir)
end

% channels = {'hoechst','bodipy','pparg','cebpb'};
name1 = 'Hoechst';
name2 = 'ddx21';
% name3 = 'pparg-mchy';
% name4 = 'pparg-mchy';

nucr = 6;
debrisarea = 100;
boulderarea = 1000;
ringcalc = 0;

shot = [num2str(row),'_',num2str(col),'_',num2str(site)];

raw1 = double(imread([rawdir,shot,'_',name1,'_',num2str(day),'.tif'])); 
raw2 = double(imread([rawdir,shot,'_',name2,'_',num2str(day),'.tif'])); 
% raw3 = double(imread([rawdir,shot,'_',name3,'_',num2str(day),'.tif'])); 
% raw4 = double(imread([rawdir,shot,'_',name4,'_',num2str(day),'.tif'])); 


% raw1 = raw1-bgcmos;
raw1 = imtophat(raw1,strel('disk',25));

% raw2 = raw2-bgcmos;
raw2 = imtophat(raw2,strel('disk',25));

% raw3 = raw3-bgcmos;
% raw3 = imtophat(raw3,strel('disk',25));
% 
% raw4 = raw3; 

raw1gray = mat2gray(raw1);
thresh = minerrthresh(raw1gray);
nuc_mask = raw1gray>thresh;
% sigmask = sigreal>minerrthresh(sigreal);

foreground = nuc_mask;

% nucmask = markershed(nucmask,2);
nuc_mask = markershed_regionalmax(nuc_mask,raw1,6);
nuc_mask = bwareaopen(nuc_mask,debrisarea);
% 
splitlength = 4;       
nuc_mask = segmentdeflections_bwboundaries(nuc_mask,splitlength,debrisarea);
[nuc_mask,~] = excludelargeandwarped(nuc_mask,boulderarea);

blurradius = 3;
real1 = raw1;
% blur1=imfilter(raw1,fspecial('disk',blurradius),'symmetric');
real2=imfilter(raw2,fspecial('disk',blurradius),'symmetric');
% blur3=imfilter(raw3,fspecial('disk',blurradius),'symmetric');
% blur4=imfilter(raw4,fspecial('disk',blurradius),'symmetric');

[nuc_label,numcells]=bwlabel(nuc_mask);

nuc_info=struct2cell(regionprops(nuc_label,real1,'Area','Centroid','MeanIntensity')');
nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';
nuc_density=squeeze(cell2mat(nuc_info(3,1,:)));
nuc_mass=nuc_density.*nuc_area;
nuc_info=regionprops(nuc_label,'PixelIdxList');
% nuc_label_extend = bwmorph(nuc_label,'thicken',3);
% nuc_info_extended = regionprops(nuc_label_extend,'PixelIdxList');
sig1 = zeros(numcells,1);
sig2 = zeros(numcells,1); %sig2tot = zeros(numcells,1);
% sig3 = zeros(numcells,1); %sig3tot = zeros(numcells,1);
% sig4 = zeros(numcells,1);
% sig4median = zeros(numcells,1); sig4mean = zeros(numcells,1);
% sig3 = zeros(numcells,1);
nanvec=ones(numcells,1)*NaN;
for cell=1:numcells
sig1(cell)=mean(real1(nuc_info(cell).PixelIdxList));

sig2(cell)=mean(real2(nuc_info(cell).PixelIdxList));
% sig2tot(cell) = mean(real2(nuc_info_extended(cell).PixelIdxList))*numel(nuc_info(cell).PixelIdxList);
% sig2mean(cell) = mean(real2(nuc_info(cell).PixelIdxList));

% sig3(cell)=median(real3(nuc_info(cell).PixelIdxList));
% sig3tot(cell) = mean(real3(nuc_info_extended(cell).PixelIdxList))*numel(nuc_info(cell).PixelIdxList);
% sig3mean(cell) = mean(real3(nuc_info(cell).PixelIdxList));
% sig4(cell)= median(real4(nuc_info(cell).PixelIdxList));
% sig4median(cell)=median(real4(nuc_info(cell).PixelIdxList));
% sig4mean(cell) = mean(real4(nuc_info(cell).PixelIdxList));
end
if ringcalc==1
        innerrad=1; outerrad=10; %10xB1|20xB2: 1/5
        ring_label=getcytoring_thicken(nuc_label_all,innerrad,outerrad,real2);
        ring_label_all = imfill(ring_label,'holes');
        ring_info=regionprops(ring_label,'PixelIdxList');   
        ring_info_all = regionprops(ring_label_all,'PixelIdxList');
        sig2ring_75th=nanvec; sig2ring_fgmedian=nanvec; sig2ring_fgtotal=nanvec;
        sig3ring_75th=nanvec; sig3ring_fgmedian=nanvec; sig3ring_fgtotal=nanvec;        
        for cell=1:numcells
            if cell>numel(ring_info)
                break;
            end
            
            ring2all=real2(ring_info(cell).PixelIdxList);
            ring2all(ring2all>prctile(ring2all,99))=[];
            sig2ring_75th(cell)=prctile(ring2all,75);
            ring2cell = real2(ring_info(cell).PixelIdxList);
            ring2cell(ring2cell<0) = [];
            sig2ring_fgtotal(cell) = sum(ring2cell);
            ring2foreground = ring2all(ring2all>0);
            
            ring3all=real3(ring_info(cell).PixelIdxList);
            ring3all(ring3all>prctile(ring3all,99))=[];
            sig3ring_75th(cell)=prctile(ring3all,75);
            ring3cell = real3(ring_info(cell).PixelIdxList);
            ring3cell(ring3cell<0) = [];
            sig3ring_fgtotal(cell) = sum(ring3cell);
            ring3foreground = ring3all(ring3all>0);
            
%             histogram(ring2all,25)
%             if numel(ring2foreground)<100
%                  ring2foreground=ring2all;
%             end
            if numel(ring2all)>50
                 sig2ring_fgmedian(cell)=nanmedian(ring2foreground);
%                  numbins=15;
%                  sig2ring_fgmode(cell)=getmode(ring2foreground,numbins);
            end
            if numel(ring3all)>50
                 sig3ring_fgmedian(cell)=nanmedian(ring3foreground);
%                  numbins=15;
%                  sig2ring_fgmode(cell)=getmode(ring2foreground,numbins);
            end            
        end
end

% goodcellflag = [ones(numgoodcells,1);zeros(numcells-numgoodcells,1)];

% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,repmat(numcells,numel(sig3),1),goodcellflag];
fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig2ring_fgtotal,sig3,sig3ring_75th,sig3ring_fgmedian,sig3ring_fgtotal,sig4,goodcellflag];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig2ring_75th,sig2ring_fgmedian,sig3,goodcellflag];
% fixdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,nuc_mass,sig1,sig2,sig3,sig4,goodcellflag];
save([datadir,'fixdata_','Plate_',num2str(day),'_',shot,'.mat'],'fixdata');
disp([shot,'_Plate_',num2str(day)])

% save([datadir,'fixdata_Day4_',shot,'.mat'],'fixdata');
% disp(shot)
toc
% end

