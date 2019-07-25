clear all; clc; close all; fclose all; direc = pwd;
%% ------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% author : JAWAD FAYAZ (email: jfayaz@uci.edu)
%  visit: (https://jfayaz.github.io)

% ------------------------------ Instructions ------------------------------------- 
% This code develops the RotD50 Sa and RotD100 Sa Spectra of the Bi-Directional 
% Ground Motion records provided in the 'GM' folder which must be in current folder. 
% 
% MAKE SURE TO HAVE OPENSEES TERMINAL IN THIS DIRECTORY BEFORE RUNNING THIS.
% MAKE SURE NOT TO HAVE ANY SPACE IN THE ENTIRE PATH TO THIS CODE (if so OpenSees
% will give an error "WARNING - PathSeries::PathSeries() - could not open file")
%
% The two directions of the ground motion record must be named as 'GM1i' and 'GM2i',
% where 'i' is the ground motion number which goes from 1 to 'n', 'n' being the total
% number of ground motions for which the Spectra needs to be generated. The extension
% of the files must be '.AT2'
% 
% For example: If the Spectra of two ground motion records are required, 4 files with
% the following names must be provided in the given 'GM' folder:
%     'GM11.AT2' - Ground Motion 1 in direction 1 (direction 1 can be either one of the bi-directional GM as we are rotating the ground motions it does not matter) 
%     'GM21.AT2' - Ground Motion 1 in direction 2 (direction 2 is the other direction of the bi-directional GM)
%     'GM12.AT2' - Ground Motion 2 in direction 1 (direction 1 can be either one of the bi-directional GM as we are rotating the ground motions it does not matter)  
%     'GM22.AT2' - Ground Motion 2 in direction 2 (direction 2 is the other direction of the bi-directional GM)
% 
% The Ground Motion file must be a vector file with 4 header lines.The first 3 lines can have
% any content, however, the 4th header line must be written exactly as per the following example:
%     'NPTS=  15864, DT= 0.0050'
%     
%
% INPUT:
% This codes provides the option to have 3 different regions of developing the Spectra of ground motions with different period intervals (discretizations)
% The following inputs within the code are required:
% 
%     'Int_T_Reg_1'        --> Period Interval for the first region of the Spectrum 
%     'End_T_Reg_1'        --> Last Period of the first region of the Spectrum (where to end the first region)
%     'Int_T_Reg_2'        --> Period Interval for the second region of the Spectrum 
%     'End_T_Reg_2'        --> Last Period of the second region of the Spectrum (where to end the second region)
%     'Int_T_Reg_3'        --> Period Interval for the third region of the Spectrum 
%     'End_T_Reg_3'        --> Last Period of the third region of the Spectrum (where to end the third region)
% 
%     'Plot_Spectra'       --> whether to plot the generated Spectra of the ground motions (options: 'Yes', 'No')    
% 
% 
% OUTPUT:
% The output will be provided in a saperate 'GMi_Spectra.txt' file for each ground motion record, where 'i' denotes the number of ground motion in the same of
% provided 'GM1i.AT2' and 'GM2i.AT2' files. The output files will be generated in a saperate folder 'Spectra' which will be created in the current folder
% The 'GMi_Spectra.txt' file will consist of space-saperated file with:
%     'Periods (secs)' 'RotD50 Sa (g)' 'RotD100 Sa (g)' 
% 
% Note: During the process another folder 'Results' will be created and
% deleted automatically at the end.
%%-------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ======================== USER INPUTS =============================== %%

% For periods 0 to 'End_T_Reg_1' in an interval of 'Int_T_Reg_1' 
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_1' = 0.01
Int_T_Reg_1   = 0.02;
End_T_Reg_1   = 1;

% For periods ['End_T_Reg_1'+'Int_T_Reg_2'] to 'End_T_Reg_2' in an interval of 'Int_T_Reg_2'
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_2' = 0.01
Int_T_Reg_2   = 0.1;
End_T_Reg_2   = 2;

% For periods ['End_T_Reg_2'+'Int_T_Reg_3'] to 'End_T_Reg_3' in an interval of 'Int_T_Reg_3'
% (Note interval must not be less than 2 decimal points, which means min 'Int_T_Reg_3' = 0.01
Int_T_Reg_3   = 0.25;
End_T_Reg_3   = 5;

% Plot Spectra  (options: 'Yes' or 'No')
Plot_Spectra  = 'Yes';


%%%%%%================= END OF USER INPUT ========================%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% =============== Running OpenSees ================
T1_Reg = [Int_T_Reg_1:Int_T_Reg_1:End_T_Reg_1];
T2_Reg = [End_T_Reg_1+Int_T_Reg_2:Int_T_Reg_2:End_T_Reg_2];
T3_Reg = [End_T_Reg_2+Int_T_Reg_3:Int_T_Reg_3:End_T_Reg_3];
Periods= [T1_Reg T2_Reg T3_Reg]';                                       

No_of_GMs = length(dir([direc,'\GMs\*.AT2']))/2;

txt = fileread('Generate_Spectra.tcl');
data = strsplit(txt,'\n')';
data{23} = ['set No_of_GMs ',round(num2str(No_of_GMs))];
data{27} = ['for {set k 1} {$k <= ',round(num2str(length(Periods))),'} {incr k 1} { '];
data{29} = ['for {set i ',num2str(round(Int_T_Reg_1*100)),'} {$i <= ',num2str(round(End_T_Reg_1 *100)),'} {incr i ',num2str(round(Int_T_Reg_1*100)),'} {'];
data{33} = ['for {set i ',num2str(round((End_T_Reg_1+Int_T_Reg_2)*100)),'} {$i <= ',num2str(round(End_T_Reg_2 *100)),'} {incr i ',num2str(round(Int_T_Reg_2*100)),'} {'];
data{37} = ['for {set i ',num2str(round((End_T_Reg_2+Int_T_Reg_3)*100)),'} {$i <= ',num2str(round(End_T_Reg_3 *100)),'} {incr i ',num2str(round(Int_T_Reg_3*100)),'} {'];
fid = fopen('Generate_Spectra.tcl','w');
fprintf(fid,'%s\n',data{:});
fclose(fid);

fprintf('\n Starting OpenSees Runs...\n')
!OpenSees.exe Generate_Spectra.tcl
fprintf('\n OpenSees Finished!\n')


%% ============= Generating Spectra using OpenSees Results ================
mkdir Spectra                     

for GM = 1:No_of_GMs
    i = 0;
    cd (['Results\RESULTS',num2str(GM)])
    for j = 1:length(Periods)
        T = Periods(j);
        i = i +1;
        if T == round(T)
            file = ['Node',num2str(GM),'_',num2str(T),'.0.out'];
        else
            file = ['Node',num2str(GM),'_',num2str(T),'.out'];
        end
        data = dlmread(file);
        x_y = [data(:,2),data(:,3)];
        omega = 2*pi/T;

        for theta = 0:1:180
            Rot_Matrix (1,1) = cosd(theta);
            Rot_Matrix (1,2) = -sind(theta);
            Rot_Matrix (2,1) = sind(theta);
            Rot_Matrix (2,2) = cosd(theta);
            Rot_Disp{theta+1,1} = x_y*Rot_Matrix;    
            clearvars Rot_Matrix
        end
        
        for l = 1:length(Rot_Disp)
            disp = Rot_Disp{l,1};
            max_disp(l,1) = max(abs(disp(:,1)));
        end
        
        Acc = max_disp.*((omega^2)/386.1);
        
        ROTD100(i,1) = max(Acc); 
        ROTD50(i,1) = median(Acc); 
       clearvars data Disp Acc disp max_disp Rot_Disp
    end
    cd (direc);
    
    if strcmpi(Plot_Spectra,'Yes') == 1
        figure(1)
        plot(Periods,ROTD100,'LineWidth',2)
        xlabel ('Period (sec)','fontsize',16,'fontWeight','bold');
        ylabel ('RotD100 Sa (g)','fontsize',16,'fontWeight','bold');
        title ('RotD100 Spectra','fontsize',16,'fontWeight','bold');
        ylim([0 ceil(max(ROTD100))])
        xlim([0 ceil(max(Periods))])
        set(gca,'fontsize',14,'FontName', 'Times New Roman','LineWidth', 1.5)
        grid on; box off
        hold on

        figure(2)
        plot(Periods,ROTD50,'LineWidth',2)
        xlabel ('Period (sec)','fontsize',16,'fontWeight','bold');
        ylabel ('RotD50 Sa (g)','fontsize',16,'fontWeight','bold');
        title ('RotD50 Spectra','fontsize',16,'fontWeight','bold');
        ylim([0 ceil(max(ROTD50))])
        xlim([0 ceil(max(Periods))])
        set(gca,'fontsize',14,'FontName', 'Times New Roman','LineWidth', 1.5)
        grid on; box off
        hold on
    end
    
    cd ([direc,'\Spectra']);
    fid = fopen (['GM',num2str(GM),'_Spectra.txt'],'w');
    fprintf(fid,'Period(s) RotD50Sa(g) RotD100Sa(g)\n');
    fprintf (fid,'%.2f %.6f %.6f\n', [Periods,ROTD50,ROTD100]');
    fclose(fid); 
    cd (direc);
    
    fprintf(' Spectrum generated for Ground Motion %d\n',GM)
    clearvars -except direc i j Plot_Spectra No_of_GMs trial all_direc Periods
end
cd (direc);
rmdir Results s