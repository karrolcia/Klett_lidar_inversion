clear all
filename1 = 'fire_27_merged_with_mask_v7.dat';
string = filename1;
outfile1 =  [ string(1:end-4) '.mat'];
fid1 = fopen(filename1);
% read the file headers, find nx (one value)
nx = fscanf(fid1, '"nx="\t%d\n', 1); % number of observations - timesteps
nz = fscanf(fid1, '"nz="\t%d\n',1); % number of levels - heights
a = textscan(fid1, '%f', 1, 'HeaderLines', 8);

for i = 1:nx
    DATA_sim {i,:} = textscan(fid1, '', nz, 'CommentStyle', '"','CollectOutput',1 );
end 

for i = 1:nx
    
    height (:, i) = DATA_sim{i,1}{1,1}(2:nz,1);   % in km
    lid_total (:,i) = DATA_sim{i,1}{1,1}(2:nz,2); % in [1/m/sr]
    lin_depol (:,i) = DATA_sim{i,1}{1,1}(2:nz,3); % 
    Ze (:,i) = DATA_sim{i,1}{1,1}(2:nz,4); % in [mm^6/m^3]
    Vd (:,i) = DATA_sim{i,1}{1,1}(2:nz,5); % in [m/s]
    T (:,i) = DATA_sim{i,1}{1,1}(2:nz,6); % in K
    P (:,i) = DATA_sim{i,1}{1,1}(2:nz,7); % in mb
    WV (:,i) = DATA_sim{i,1}{1,1}(2:nz,8); % in [g/m^3]
    MASKb (:,i) = DATA_sim{i,1}{1,1}(2:nz,9); 
    LWP (:,i) = DATA_sim{i,1}{1,1}(2:nz,10); % in cm
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


filename2 = 'true_fire_27_merged_with_mask_v7.dat';
string2 = filename2;
outfile2 =  [ string2(1:end-4) '.mat'];
fid2 = fopen(filename2);
% read the file headers, find nx (one value)
% nx = fscanf(fid, '"nx="\t%d\n', 1);
% nz = fscanf(fid, '"nz="\t%d\n',1);
a = textscan(fid2, '%f', 1, 'HeaderLines', 9);

for i = 1:nx
    DATA_real {i,:} = textscan(fid2, '', nz, 'CommentStyle', '"','CollectOutput',1 );
end 

for i = 1:nx
    
    height_real (:, i) = DATA_real{i,1}{1,1}(2:nz,1);   % in km
    No (:,i) = DATA_real{i,1}{1,1}(2:nz,2); % in [1/cm3]
    Ext (:,i) = DATA_real{i,1}{1,1}(2:nz,3); % in [1/m]
    LWC (:,i) = DATA_real{i,1}{1,1}(2:nz,4); % in [g/m3] 
   
end

% Transposing vriables to nt x nzg
height_real = height_real' ;
No = No' ;
Ext = Ext';
LWC = LWC';

filename3 = 'fire_27_merged_with_mask_v7_ss.dat';
string3 = filename3;
outfile3 =  [ string3(1:end-4) '.mat'];
fid3 = fopen(filename3);
% read the file headers, find nx (one value)
% nx = fscanf(fid, '"nx="\t%d\n', 1);
% nz = fscanf(fid, '"nz="\t%d\n',1);
a = textscan(fid3, '%f', 1, 'HeaderLines', 9);

for i = 1:nx
    DATA_sim_ss {i,:} = textscan(fid3, '', nz, 'CommentStyle', '"','CollectOutput',1 );
end 

for i = 1:nx
    height_ss (:, i) = DATA_sim_ss{i,1}{1,1}(2:nz,1);   % in km
    lid_total_ss (:,i) = DATA_sim_ss{i,1}{1,1}(2:nz,2); % in [1/m/sr]
    lin_depol_ss (:,i) = DATA_sim_ss{i,1}{1,1}(2:nz,3); % 

   
end
height_ss = height_ss';
lid_total_ss = lid_total_ss';
lin_depol_ss = lin_depol_ss';

% filename4 = 'true_fire_27_merged_with_mask_v7_ss.dat';
% string4 = filename4;
% outfile4 =  [ string4(1:end-4) '.mat'];
% fid2 = fopen(filename4);
% % read the file headers, find nx (one value)
% % nx = fscanf(fid, '"nx="\t%d\n', 1);
% % nz = fscanf(fid, '"nz="\t%d\n',1);
% a = textscan(fid2, '%f', 1, 'HeaderLines', 9);
% 
% for i = 1:nx
%     DATA_real {i,:} = textscan(fid2, '', nz, 'CommentStyle', '"','CollectOutput',1 );
% end 
% 
% for i = 1:nx
%     
%     height_real_ss (:, i) = DATA_real{i,1}{1,1}(2:nz,1);   % in km
%     No_ss (:,i) = DATA_real{i,1}{1,1}(2:nz,2); % in [1/cm3]
%     Ext_ss (:,i) = DATA_real{i,1}{1,1}(2:nz,3); % in [1/m]
%     LWC_ss (:,i) = DATA_real{i,1}{1,1}(2:nz,4); % in [g/m3] 
%    
% end
% 
% % Transposing vriables to nt x nzg
% % height_real_ss = height_real' ;
% % No_ss = No' ;
% % Ext_ss = Ext';
% % LWC_ss = LWC';
% 
% % Transposing vriables to nt x nzg
% 
% 
save (outfile1, 'height' , 'lid_total', 'lin_depol', 'Ze', 'Vd', 'T' , 'P', ...
    'WV','LWP' ,'nx','nz')
save (outfile2, 'height_real', 'No', 'Ext' , 'LWC')
save (outfile3, 'height_ss' , 'lid_total_ss', 'lin_depol_ss')
% save (outfile4, 'height_real_ss', 'No_ss', 'Ext_ss' , 'LWC_ss')
