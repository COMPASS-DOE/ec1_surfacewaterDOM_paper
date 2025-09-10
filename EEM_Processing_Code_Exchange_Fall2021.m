%% Section 1: Important things to change prior to running the code

mac = false

%Put any important Notes about the runs that are relevant here
special_notes='Fluorescence data was blank subtracted, inner filter corrected, and raman normalized with instrument software. Exported from MatLab on June 16, 2022.'

%Set Parameters for manual plotting
nrows = 1;
ncols = 1;
clim_min = 0; %minimum fluorescence intensity graphed
clim_max = 0.45 %maximum fluorescence intensity graphed;  use max(eem.X(:),[],'omitnan') to normalize to the maximum peak intesity within the whole dataset
manualscale = false; %Set as true if you want to use the minimum/maximum scales set above. This will apply to all eems. Set to false if you want to autoscale for each individual eem 
autosave = true; %True will autosave all eems; false will no save and will display images individually
figuretitle = true %Will put a title at the top of the EEM based on the Sample Description tab in the sample log

%Is DOC concentration available in the sample log? Put true if yes, put
%false if no
DOCconcentration = true
DOC_normalized_figures = true

%Demo - put a copy of the demo dataset in a writeable folder, e.g.
demopath='C:\Users\roeb782\OneDrive - PNNL\Desktop\EEMs\Processed Online\2021 Runs\EXCHANGE_Fall2021\PARAFAC'; %edit for correct file destination


studycode = 'Exchange_2021Oct14';
 


%% Section 2: Reading in EEMs and Absorbance Files
cd([demopath '/2. Instrument EEMs'])
filetype=3;ext = 'dat';RangeIn='A4..GJ253';headers=[0 1];display_opt=0;outdat=0;
[X,Emmat,Exmat,filelist_eem,outdata]=readineems(filetype,ext,RangeIn,headers,display_opt,outdat);
Ex=Exmat(1,:); %Since all files have the same excitation wavelengths
Em=Emmat(:,1); %Since all files have the same emission wavelengths


cd([demopath '/1. Absorbance'])
filetype='Abs';ext = 'dat';RangeIn='A2..B192';display_opt=0;outdat=0;
[S_abs,W_abs,wave_abs,filelist_abs]=readinscans(filetype,ext,RangeIn,display_opt,outdat);

%% Section 3: Reading in Sample Log

% Sample Log
cd(demopath)
[LogNUM,LogTXT]=xlsread('SampleLog.xlsx','log');

%Get filenames and text data from LogTXT
Log_Date=LogTXT(:,2);
Log_SampleDescription=LogTXT(:,3);
Log_SampleID=LogTXT(:,4);
log_IntegrationTime=LogTXT(:,6);
Log_EEMfile=LogTXT(:,7);
Log_ABSfile=LogTXT(:,8);

%Get numeric data from LogNUM
RepNo=LogNUM(:,5);
pathlength=LogNUM(:,9);
df=LogNUM(:,10);
if DOCconcentration
DOCconc=LogNUM(:,11);
end

%% Section 4: Match EEM and Absorbance Files with information from Sample Log

%Pair the EEMs with the other datasets
Pair_EEM_log=[Log_EEMfile Log_EEMfile]; %used for all numeric information in the log
Pair_EEM_abs=[Log_EEMfile Log_ABSfile];

% Obtain matching datasets - from loaded datasets
Sabs=matchsamples(filelist_eem,filelist_abs,Pair_EEM_abs,X,S_abs);% ABS scans that match filelist_eem

% numbers only                                          
replicates=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,RepNo);      % RepNo matching filelist_eem
dilfac=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,df);             % df matching filelist_eem
Abspath=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,pathlength);
if DOCconcentration
DOC=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,DOCconc);
end

% text only - need to remove headers, e.g. Log_Site(2:end,:);  
sample_description=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_SampleDescription(2:end,:));
sampleID=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_SampleID(2:end,:));     % sites matching filelist_eem
dates=matchsamples(filelist_eem,Log_EEMfile(2:end,:),Pair_EEM_log,X,Log_Date(2:end,:));         % dates matching filelist_eem


%% Section 5: Correct Absorbance and EEM files for pathlength, dilutions, and DOC normalizations

%Correct Absorbance Data for Pathlength
Sabs_pl=Sabs./Abspath

%Create Absorbance 2D matrix
A=[wave_abs;Sabs_pl]; 

% Eliminate data below or above the wavelength range of Ex and Em correction files
Em_in=Em(Em<=600); 
X_in=X(:,Em<=600,:); 

% Optional - correct for dilution if necessary
% If samples were diluted before measuring EEMs and Abs scans, divide the
% corrected EEMs by the sorted dilution factors (dilfac NOT df).
%Dilution Correct EEMs and Absorbance
XcRUn_df=X_in./dilfac %Dilution Corrects EEMs
Sabs_pl_df=Sabs_pl./dilfac %Dilution Corrects Absorbance 

%DOC normalization to 1 mg/L
if DOCconcentration
XcRUn_df_doc=XcRUn_df./DOC %Normalize EEMs to 1 mg/L DOC
Sabs_pl_df_doc=Sabs_pl_df./DOC %Normalize Absorbance to 1 mg/L DOC
else
end


%% Section 6: Building Final Sample Set
%Build dataset with Raw Data (e.g. not DOC normalized)
mydata_raw=assembledataset(XcRUn_df,Ex,Em_in,'RSU','filelist',filelist_eem,'sampleID',sampleID,'sample_description',sample_description,'rep',replicates,[4]) %Build new 3D Matrix

%Build dataset with DOC normalized data
if DOCconcentration
mydata_doc=assembledataset(XcRUn_df_doc,Ex,Em_in,'RSU','filelist',filelist_eem,'sampleID',sampleID,'sample_description',sample_description,'rep',replicates,[4]) %Build new 3D Matrix
else
end

if DOCconcentration
    if DOC_normalized_figures
        SubData=subdataset(mydata_doc,[],mydata_doc.Em>550,mydata_doc.Ex<260);
        SubData=subdataset(SubData,[],SubData.Em<260,SubData.Ex>450)
    else
        SubData=subdataset(mydata_raw,[],mydata_raw.Em>550,mydata_raw.Ex<260);
        SubData=subdataset(SubData,[],SubData.Em<260,SubData.Ex>450)
    end
else
    SubData=subdataset(mydata_raw,[],mydata_raw.Em>550,mydata_raw.Ex<260);
    SubData=subdataset(SubData,[],SubData.Em<260,SubData.Ex>450)
end

Xs=smootheem(SubData,[30 10],[3 45],[50 30],[50 10],[1 1 0 0],[20],3500,''); %Removes 1st and 2nd order Raman Scatter. Attempts to smooth 1st order. Note Rayleigh Scatter removed during instrument processing. 

eemview(Xs,'X',[1 1],30,[],{'sampleID'},[],[],'colorbar')

%% Section 7: Saving Final Figures and Exporting Spectral Indices File
cd([demopath '/5. Figures'])
manualplots


cd([demopath '/6. Spectral Indices'])
SpectralIndicesExport_Updated_July2022(1)

%% Section 8: Export Final EEMs Data

cd([demopath '/3. Corrected_EEMs'])
Xc = permute(XcRUn_df,[2 3 1]); %Rearranges 3D matrix containing dilution and raman normalized data in a way so it can be saved to individual files. 
Xc2=flip(Xc,2); %flips data along the 2nd diminsion, puts data in the same format as the Aqualog input files
Exflip=fliplr(Ex)% flips excitation wavelengths to be consistent with Aqualog in put files
Ex_cell=num2cell(Exflip)
Ex_new=["Wavelength" Ex_cell]% adds a 0 to excitation wavelengths in the front so that dimenstions of data and wavelengths being concatenated below are all consistent. Required to do or an error occurs. 

%This loop takes individual EEM data from the dilution corrected and raman normalized matrix, pairs with correct Excitation and Emmison wavelengths,
%and exports individual corrected EEM files as a csv. 
for i=1:size(Xc2,3)
    filename=char(sampleID(i));
    filename=[filename,'_DilCorr_IFE_RamNorm.dat'];
    X3=Xc2(:,:,i);
    C = [Em_in X3];
    D = [Ex_new;C];
    writematrix(D,filename,'Delimiter','tab') 
end

%% Section 9: Export Final EEMs Data
cd([demopath '/4. Dil_Corrected_Absorbance'])
Sabs_pl_df_t = Sabs_pl_df'; %Transpose dilution corrected absorbance data
wave_abs_t=wave_abs'

for i=1:size(Sabs_pl_df_t,2)
    filename=char(sampleID(i));
    filename=[filename,'_DilCorr_Abs.dat'];
    Absorbance = Sabs_pl_df_t(:,i);% Extract this one column into its own variable.
    E=[wave_abs_t Absorbance]
    F=["Wavelength (nm)" "Dilution Corrected Absorbance"];
    G=[F;E];
    writematrix(G,filename,'Delimiter','tab')
end
  

