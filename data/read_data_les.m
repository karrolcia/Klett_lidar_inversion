clear all
filename1 = 'les_1-merged_with_driz.dat';
string = filename1;
outfile1 =  [ string(1:end-4) '.mat'];
fid1 = fopen(filename1);
% read the file headers, find nx (one value)
 nx = fscanf(fid1, '"nx="\t%d\n', 1);
 nz = fscanf(fid1, '"nz="\t%d\n',1);
a = textscan(fid1, '%f', 1, 'HeaderLines', 6);


for i = 1:nx
    DATA_sim_les {i,:} = textscan(fid1, '', nz, 'CommentStyle', '"','CollectOutput',1 );
end 

for i = 1:nx
    
    height (:, i) = DATA_sim_les{i,1}{1,1}(1:160,1);   % in km
    lid_total (:,i) = DATA_sim_les{i,1}{1,1}(1:160,2); % in [1/m/sr]
    lin_depol (:,i) = DATA_sim_les{i,1}{1,1}(1:160,3); % 
    Ze (:,i) = DATA_sim_les{i,1}{1,1}(1:160,4); % in [mm^6/m^3]
    Vd (:,i) = DATA_sim_les{i,1}{1,1}(1:160,5); % in [m/s]
    T (:,i) = DATA_sim_les{i,1}{1,1}(1:160,6); % in K
    P (:,i) = DATA_sim_les{i,1}{1,1}(1:160,7); % in mb
    WV (:,i) = DATA_sim_les{i,1}{1,1}(1:160,8); % in [g/m^3]
    MASKb (:,i) = DATA_sim_les{i,1}{1,1}(1:160,9); 
    LWP (:,i) = DATA_sim_les{i,1}{1,1}(1:160,10); % in cm
    total (:,i) = DATA_sim_les{i,1}{1,1}(1:160,25); % in [1/m/sr]
    Ray (:,i) = DATA_sim_les{i,1}{1,1}(1:160,26); % in [1/m/sr]
    total2 (:,i) = DATA_sim_les{i,1}{1,1}(1:160,27); % in [1/m/sr]
    Ray2 (:,i) = DATA_sim_les{i,1}{1,1}(1:160,28); % in [1/m/sr]
end

% Transposing vriables to nt x nz
height = height';
lid_total = lid_total' ;
lin_depol = lin_depol' ;
Ze = Ze' ;
Vd = Vd' ;
T = T' ;
P = P';
WV = WV';
MASKb= MASKb' ;
LWP = LWP' ;
total = total';
total2 = total2' ;
Ray = Ray';
Ray2 = Ray2' ;



filename2 = 'true_les_1-merged_with_driz.dat';
string2 = filename2;
outfile2 =  [ string2(1:end-4) '.mat'];
fid2 = fopen(filename2);
% read the file headers, find nx (one value)
% nx = fscanf(fid, '"nx="\t%d\n', 1);
% nz = fscanf(fid, '"nz="\t%d\n',1);
a = textscan(fid2, '%f', 1, 'HeaderLines', 7);

for i = 1:nx
    DATA_real_les {i,:} = textscan(fid2, '', nz, 'CommentStyle', '"','CollectOutput',1 );
end 

for i = 1:nx
    
    height_real (:, i) = DATA_real_les{i,1}{1,1}(1:160,1);   % in km
    No (:,i) = DATA_real_les{i,1}{1,1}(1:160,2); % in [1/cm3]
    Ext (:,i) = DATA_real_les{i,1}{1,1}(1:160,3); % in [1/m]
    LWC (:,i) = DATA_real_les{i,1}{1,1}(1:160,4); % in [g/m3] 
   
end

% Transposing vriables to nt x nz
height_real = height_real' ;
No = No' ;
Ext = Ext';
LWC = LWC';

save (outfile1, 'height' , 'lid_total', 'lin_depol', 'Ze', 'Vd', 'T' , 'P', ...
    'WV','LWP','Ray','nx','nz')
save (outfile2, 'height', 'No', 'Ext' , 'LWC')