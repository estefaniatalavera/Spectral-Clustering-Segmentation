%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%
%
% Read Spectral Clustering Results & Plot 

%close all
addpath('../Results/SpectralClustering/ResultsFM&JI_511')
loadParameters;

%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.

for ind_dataset=1:length(Dataset_ops)
    
Dataset=Dataset_ops{ind_dataset};
if strcmp(Dataset,'Narrative')==1
    
    folders_indx={'Petia1','Petia2','Mariella','Estefania1','Estefania2'};

elseif strcmp(Dataset,'SenseCam')==1
    
    folders_indx={'Day1','Day2','Day3','Day4','Day6'};

elseif strcmp(Dataset,'All')==1
    
	folders_indx={'Petia1','Petia2','Mariella','Estefania1','Estefania2','Day1','Day2','Day3','Day4','Day6'};
 
end

%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%
%Plot parameters
colorlist={'r','b','m','g','c','k','y'};
nColours = 3; % aqui cambias cuantos colores quieres
c = colormap(jet);  % aqui el colormap
close(gcf); 
colours = c(round(linspace(1,size(c,1),nColours)),:);



%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%
% FMEASURE
for i_matrix=1:length(sim_matrix)
    
    figure,
    sim_matrix_act=sim_matrix{i_matrix};
        if strcmp(sim_matrix_act,'NN')==1,
            Spectral_Param=NN;
        elseif strcmp(sim_matrix_act,'Sigma')==1,
            Spectral_Param=Sig;            
        elseif strcmp(sim_matrix_act,'Epsilon')==1,
            Spectral_Param=Eps;            
        end    
    
    for Type=1:1:3
    
        disp(['Metodo ' sim_matrix_act ' & ' num2str(Type)])
    
        %Mean de todos los folders
        vect=[];        array=[];
        for i_fold=1:length(folders_indx)
            folder=folders_indx{i_fold};
            disp(['Metodo ' sim_matrix_act ' & folder ' folder ])

            %Load Spectral Clustering Results
            load(['Res_Spec_' folder '_' sim_matrix_act '_ParamVal_' num2str(Spectral_Param) '_Type_' num2str(Type) '.mat']);

            for row=1:size(Results,2), 
                array(row,1)=Results{row}.RPAF_Spectral.fMeasure;
            end

            vect(:,i_fold)=array';
        end%folder 

        %Mean por col (vect= [filas=cutvalu,col=folders])     
        vect_mean=mean(vect,2); 
        e=var(vect');

        %Plot mean FMeasure
        errorbar(k_values,vect_mean',e','Color',colorlist{Type},'Linewidth',3); 
        xlim([min(k_values),max(k_values)]); ylim([0,1]);
        hold on
    
    end %end type

    title(['fMeasure ' Dataset ' ' sim_matrix_act], 'FontSize', 20); 
    xlabel('cut Value', 'FontSize', 20);
    ylabel('f-Measure', 'FontSize', 20);
    legend('Type 1','Type 2','Type 3');
    set(gca,'FontSize',20); grid on ,
    disp('Done JIMean');
    
end%end Fm



%% -.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
%
% Jaccard Index
for i_matrix=1:length(sim_matrix)
    
    figure,
    sim_matrix_act=sim_matrix{i_matrix};
    
        if strcmp(sim_matrix_act,'NN')==1,
            Spectral_Param=NN;
        elseif strcmp(sim_matrix_act,'Sigma')==1,
            Spectral_Param=Sig;            
        elseif strcmp(sim_matrix_act,'Epsilon')==1,
            Spectral_Param=Eps;            
        end
    
    for Type=1:1:3
    
        disp(['Metodo ' sim_matrix_act ' & ' num2str(Type)])
    
        %Mean de todos los folders
        vect=[];        array=[];
        for i_fold=1:length(folders_indx)
            folder=folders_indx{i_fold};
            disp(['Metodo ' sim_matrix_act ' & folder ' folder ])

            %Load Spectral Clustering Results
            load(['Res_Spec_' folder '_' sim_matrix_act '_ParamVal_' num2str(Spectral_Param) '_Type_' num2str(Type) '.mat']);

            for row=1:size(Results,2), 
                array(row,1)=Results{row}.RPAF_Spectral.JaccardIndex;
            end

            vect(:,i_fold)=array';
        end%folder 

        %Mean por col (vect= [filas=cutvalu,col=folders])     
        vect_mean=mean(vect,2); 
        e=var(vect');

        %Plot mean FMeasure
        errorbar(k_values,vect_mean',e','Color',colorlist{Type},'Linewidth',3); 
        xlim([min(k_values),max(k_values)]); ylim([0,1]);
        hold on
    
    end %end type

    title(['Jaccard Index ' Dataset ' ' sim_matrix_act], 'FontSize', 20); 
    xlabel('cut Value', 'FontSize', 20);
    ylabel('Jaccard Index', 'FontSize', 20);
    legend('Type 1','Type 2','Type 3');
    set(gca,'FontSize',20); grid on ,
    disp('Done JIMean');
    
end%end Jaccard

end % dataset