

%% Loads parameters
loadParameters;

%% Folders
for i_fold=1:length(folders)
    
    %% Build paths for images, excel, features and results
    folder=folders{i_fold};
    fichero=([directorio_im '/' camera{i_fold} '/imageSets/' folder]);
    path_excel = [directorio_im '/' camera{i_fold} '/GT/GT_' folder '.xls'];
    path_features = [directorio_im '/' camera{i_fold} '/CNNfeatures/CNNfeatures_' folder '.mat'];
    path_features_PCA = [directorio_im '/' camera{i_fold} '/CNNfeatures/CNNfeaturesPCA_' folder '.mat'];
    root_results = [directorio_results '/' folder];
    mkdir(root_results);
    

    %% Images
    files_aux=dir([fichero '/*' formats{i_fold}]);
    count = 1;
    for n_files = 1:length(files_aux)
        if(files_aux(n_files).name(1) ~= '.')
            files(count) = files_aux(n_files);
            count = count+1;
        end
    end
    Nframes=length(files);
    

    %% Excel
    [clust_man,clustersIdGT,cl_limGT, ~]=analizarExcel_Narrative(path_excel, files);
    delim=cl_limGT';
    if delim(1) == 1, delim=delim(2:end); end
    for i=1:length(clust_man)
         [a,b]=find(clustersIdGT==i);
         clust_manId{i,1}=b;
    end     
      

    %% Features
	load(path_features);
    
    % Features
    [features_norm] = signedRootNormalization(features);

    %PCA FEATURES
    if(exist(path_features_PCA) > 0)
        load(path_features_PCA);
    else
        [ featuresPCA, ~, ~ ] = applyPCA( features_norm, paramsPCA ) ; 
        save(path_features_PCA, 'featuresPCA');
    end    
    
    %% SpectralClust
    if(paramsPCA.usePCA_Spect)
        features_Sp = featuresPCA;
    else
        features_Sp = features;
    end
    
    for matrix_indx=1:length(sim_matrix)
        SimM=sim_matrix{matrix_indx};
        
        if strcmp(SimM,'NN')==1,
            sigmaNN=0.5; Type_NN=1;
            Spectral_Param=NN;
        	W = SimGraph_NearestNeighbors(features_Sp', NN, Type_NN, sigmaNN);

        elseif strcmp(SimM,'Sigma')==1,
            Spectral_Param=Sig;
            W = SimGraph_Full(features_Sp', Sig);
            
        elseif strcmp(SimM,'Epsilon')==1,
            Spectral_Param=Eps;
            W = SimGraph_Epsilon(features_Sp', Eps);    
            
        end
        
        %%SpectralClust Type
        for Type=1:1:3
            Results={};
            %%Kvalue
            for k_indx=1:length(k_values)
                k=k_values(k_indx);
                
                disp([folder ' Clusters Id - ' SimM ' k=' num2str(k) ' Type=' num2str(Type) ' & Val=' num2str(Spectral_Param)])
                [C, L, U] = SpectralClustering(W, k, Type);
                
                clustersId=[];
                for num_k=1:k
                   vect_pos=find(C(:,num_k)==1)'; 
                   clustersId(1,vect_pos)=num_k; 
                end

                %% AFTER IDs EXTRACTION
                % we will have N Id equal to the range of sigma defined
                %for num_ids=1:length(Spectral_Param)
                index=1;
                bound=[];
                for pos=1:length(clustersId)-1
                    if (clustersId(pos)~=clustersId(pos+1))>0
                        bound(index)=pos;
                        index=index+1;
                    end
                end
                if (exist('bound','var')==0)
                    bound=0;
                end

                automatic=bound;
                if automatic(1) == 1
                    automatic=automatic(2:end);
                end
               
                % clust_man & clust_auto = array of cells     
                % LH MATRIX: Nos permite aplicar el mismo criterio que hemos
                % aplicado con FMeasuer-> separar cuando imagenes consecutivas
                % con coinciden
                LH=[];clust_autoId=[];
                for i_cl=1:max(clustersId)
                    [~,val_pos]=find(clustersId==i_cl);
                    for pos_LH=1:length(val_pos)
                        LH(val_pos(pos_LH),i_cl)=1;
                    end
                end 
                [ labels_event, ~, ~ ] = getEventsFromLH(LH);
                %Agrupamos por etiqueta
                for i_lab=1:max(labels_event)
                    [~,b]=find(labels_event==i_lab);
                    clust_autoId{i_lab,1}=b;
                end 

                % Asignar el nombre de la imagen
                clust_auto_ImagName=image_assig(clust_autoId,files);
                clust_man_ImagName=image_assig(clust_manId,files);

                [rec,prec,acc,fMeasure_Spect]=Rec_Pre_Acc_Evaluation(delim,automatic,Nframes,tol);
                [JaccardIndex_result,JaccardVar,~,~,~]=JaccardIndex(clust_man_ImagName,clust_auto_ImagName);  

                RPAF_Spectral.clustersIDs = clustersId;
                RPAF_Spectral.boundaries = bound;
                RPAF_Spectral.recall = rec;
                RPAF_Spectral.precision = prec;
                RPAF_Spectral.accuracy = acc;
                RPAF_Spectral.fMeasure = fMeasure_Spect;
                RPAF_Spectral.JaccardIndex = JaccardIndex_result;
                RPAF_Spectral.JaccardVariance = JaccardVar;               
                RPAF_Spectral.NumClusters = length(clust_auto_ImagName);

                Results{k_indx}.Similarity_Matrix = SimM;
                Results{k_indx}.Similarity_parameter = Spectral_Param;
                Results{k_indx}.Type = Type;
                Results{k_indx}.k_value = k;
                Results{k_indx}.RPAF_Spectral = RPAF_Spectral; 
 
            end%kvalue
            %% SAVE
            file_save=(['Res_Spec_' folder '_' SimM '_ParamVal_' num2str(Spectral_Param) '_Type_' num2str(Type) '.mat']);
            %save([root_results '/' file_save], 'Results');
            save(['../ResultsFM&JI50.51/' file_save], 'Results');
        end%end Type 
    end%similarity matrix
end %end folder

disp('Spectral Clustering Complete')