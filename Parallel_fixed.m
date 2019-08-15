rows = [2:7];
cols = [9:11];
sites = 1:4;
plates = 1;
manualcontrol = 0;
targetdir = 'Z:\michael\Fixed Cell Experiments\2017\20170512-nb-sirna';
missingfile = [targetdir,'Data\MissingData.mat'];

numrows=length(rows);
numcols=length(cols);
numsites=length(sites);
numplates = length(plates);
shots = numrows*numcols*numsites*numplates;


if manualcontrol==1
    load(missingfile)
    manualwells = namat;
    shots=size(manualwells,1);
end
% parpool(4)
parfor shot = 1:shots
    if manualcontrol == 1
        row = manualwells(shot,1);
        col = manualwells(shot,2);
        site = manualwells(shot,3);
        plate = manualwells(shot,4);
    else
        siteidx=mod(shot,numsites);
        if siteidx==0
            siteidx=numsites;
        end
        site=sites(siteidx);

        colidx=mod(ceil(shot/numsites),numcols);    
        if colidx==0
            colidx=numcols;
        end
        col=cols(colidx);

        rowidx=mod(ceil(shot/(numcols*numsites)),numrows);
        if rowidx==0
            rowidx=numrows;
        end
        row=rows(rowidx);

        plateidx = ceil(shot/(numcols*numsites*numrows));
        plate = plates(plateidx);
    end
    
    
 try

%     Fixed_Cell_analysis_multiple_day(row,col,site,plate);
%     Fixed_Cell_2Ring(row,col,site,plate);
    Fixed_cell_hap1(row,col,site,plate);
%     Fixed_Cell_aligned(row,col,site,plate)
%     Fixed_Cell_bodipy(row,col,site,plate);
%     Fixed_Cell_Cytosol(row,col,site,plate);
%     FISH_Only(row,col,site,plate)
%     FISH_IF(row,col,site,plate)

    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'_',num2str(plate),'\n']);  
 catch
     disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site),'_Plate_',num2str(plate)]);
 end

end
