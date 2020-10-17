close all; clear; clc;
addpath(genpath('./utils/'));
addpath(genpath('./codes/'));

result_URL = './results/';
if ~isdir(result_URL)
    mkdir(result_URL);
end

db = {'IAPRTC-12','MIRFLICKR','NUSWIDE10'};
hashmethods = {'LEMON'};
%hashmethods = {'SCM-seq','DCH','LCMFH','SCRATCH','DLFH','OCMH','OLSH','LEMON'};
loopnbits = [8 16 32 64 128];

param.top_K = 2000;

for dbi = 1     :length(db)
    db_name = db{dbi}; param.db_name = db_name;

    % load dataset
    load(['./datasets/',db_name,'.mat']);
    result_name = [result_URL 'final_' db_name '_result_test' '.mat'];

    if strcmp(db_name, 'IAPRTC-12')
        param.chunksize = 2000;
        clear V_tr V_te
        X = [I_tr; I_te]; Y = [T_tr; T_te]; L = [L_tr; L_te];
        R = randperm(size(L,1));
        queryInds = R(1:2000);
        sampleInds = R(2001:end);
        
        param.nchunks = floor(length(sampleInds)/param.chunksize);
        
        XChunk = cell(param.nchunks,1);
        YChunk = cell(param.nchunks,1);
        LChunk = cell(param.nchunks,1);
        for subi = 1:param.nchunks-1
            XChunk{subi,1} = X(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            YChunk{subi,1} = Y(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            LChunk{subi,1} = L(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
        end
        XChunk{param.nchunks,1} = X(sampleInds(param.chunksize*subi+1:end),:);
        YChunk{param.nchunks,1} = Y(sampleInds(param.chunksize*subi+1:end),:);
        LChunk{param.nchunks,1} = L(sampleInds(param.chunksize*subi+1:end),:);
        
        XTest = X(queryInds, :); YTest = Y(queryInds, :); LTest = L(queryInds, :);
        clear X Y L

    elseif strcmp(db_name, 'MIRFLICKR')
        param.chunksize = 2000;
        X = [I_tr; I_te]; Y = [T_tr; T_te]; L = [L_tr; L_te];
        R = randperm(size(L,1));
        queryInds = R(1:2000);
        sampleInds = R(2001:end);
        param.nchunks = floor(length(sampleInds)/param.chunksize);
        
        XChunk = cell(param.nchunks,1);
        YChunk = cell(param.nchunks,1);
        LChunk = cell(param.nchunks,1);
        for subi = 1:param.nchunks-1
            XChunk{subi,1} = X(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            YChunk{subi,1} = Y(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            LChunk{subi,1} = L(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
        end
        XChunk{param.nchunks,1} = X(sampleInds(param.chunksize*subi+1:end),:);
        YChunk{param.nchunks,1} = Y(sampleInds(param.chunksize*subi+1:end),:);
        LChunk{param.nchunks,1} = L(sampleInds(param.chunksize*subi+1:end),:);
        
        XTest = X(queryInds, :); YTest = Y(queryInds, :); LTest = L(queryInds, :);
        clear X Y L

    elseif strcmp(db_name, 'NUSWIDE10')
        param.chunksize = 10000;
        X = [I_tr; I_te]; Y = [T_tr; T_te]; L = [L_tr; L_te];
        R = randperm(size(L,1));
        queryInds = R(1:2000);
        sampleInds = R(2001:end);
        param.nchunks = floor(length(sampleInds)/param.chunksize);
        
        XChunk = cell(param.nchunks,1);
        YChunk = cell(param.nchunks,1);
        LChunk = cell(param.nchunks,1);
        for subi = 1:param.nchunks-1
            XChunk{subi,1} = X(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            YChunk{subi,1} = Y(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
            LChunk{subi,1} = L(sampleInds(param.chunksize*(subi-1)+1:param.chunksize*subi),:);
        end
        XChunk{param.nchunks,1} = X(sampleInds(param.chunksize*subi+1:end),:);
        YChunk{param.nchunks,1} = Y(sampleInds(param.chunksize*subi+1:end),:);
        LChunk{param.nchunks,1} = L(sampleInds(param.chunksize*subi+1:end),:);
        
        XTest = X(queryInds, :); YTest = Y(queryInds, :); LTest = L(queryInds, :);
        clear X Y L
    end
    clear I_tr I_te L_tr L_te T_tr T_te
    
    
    %% Methods
    eva_info = cell(length(hashmethods),length(loopnbits));
    
    
    for ii =1:length(loopnbits)
        fprintf('======%s: start %d bits encoding======\n\n',db_name,loopnbits(ii));
        param.nbits = loopnbits(ii);
        for jj = 1:length(hashmethods)
            
            switch(hashmethods{jj})
                case 'LEMON'
                    fprintf('......%s start...... \n\n', 'LEMON');
                    OURparam = param;
                    OURparam.alpha = 10000; OURparam.beta = 10000; OURparam.theta = 1;
                    OURparam.gamma = 0.1; OURparam.xi = 1;
                    eva_info_ = evaluate_OUR(XChunk,YChunk,LChunk,XTest,YTest,LTest,OURparam);
                case 'DCH'
                    fprintf('......%s start...... \n\n', 'DCH');
                    DCHparam = param;
                    eva_info_ = evaluate_DCH(XChunk,YChunk,LChunk,XTest,YTest,LTest,DCHparam);             
                case 'SCM-seq'
                    fprintf('......%s start...... \n\n', 'SCM-seq');
                    SCM_seqparam = param;
                    eva_info_ =evaluate_SCM_seq(XChunk,YChunk,LChunk,XTest,YTest,LTest,SCM_seqparam);
                case 'SCRATCH'
                    fprintf('......%s start...... \n\n', 'SCRATCH');
                    SCRATCHparam = param;
                    eva_info_ =evaluate_SCRATCH(XChunk,YChunk,LChunk,XTest,YTest,LTest,SCRATCHparam);
                case 'LCMFH'
                    fprintf('......%s start...... \n\n', 'LCMFH');
                    LCMFHparam = param;
                    eva_info_ =evaluate_LCMFH(XChunk,YChunk,LChunk,XTest,YTest,LTest,LCMFHparam);
                case 'DLFH'
                    fprintf('......%s start...... \n\n', 'DLFH');
                    DLFHparam = param;
                    eva_info_ =evaluate_DLFH(XChunk,YChunk,LChunk,XTest,YTest,LTest,DLFHparam);
                case 'OCMH'
                    fprintf('......%s start...... \n\n', 'OCMH');
                    OCMHparam = param;
                    eva_info_ =evaluate_OCMH(XChunk,YChunk,LChunk,XTest,YTest,LTest,OCMHparam);
                case 'OLSH'
                    fprintf('......%s start...... \n\n', 'OLSH');
                    OLSHparam = param;
                    eva_info_ =evaluate_OLSH(XChunk,YChunk,LChunk,XTest,YTest,LTest,OLSHparam);
            end
            eva_info{jj,ii} = eva_info_;
            clear eva_info_
        end
    end
    
    
    %% MAP
    for ii = 1:length(loopnbits)
        for jj = 1:length(hashmethods)
            Table_ItoT_MAP(jj,ii) = eva_info{jj,ii}{param.nchunks}.Image_VS_Text_MAP;
            Table_TtoI_MAP(jj,ii) = eva_info{jj,ii}{param.nchunks}.Text_VS_Image_MAP;
            
            for kk = 1:param.nchunks
                % MAP
                Image_VS_Text_MAP{ii}{jj,kk} = eva_info{jj,ii}{kk}.Image_VS_Text_MAP;
                Text_VS_Image_MAP{ii}{jj,kk} = eva_info{jj,ii}{kk}.Text_VS_Image_MAP;
                
                % Precision VS Recall
                Image_VS_Text_recall{ii}{jj,kk,:}    = eva_info{jj,ii}{kk}.Image_VS_Text_recall';
                Image_VS_Text_precision{ii}{jj,kk,:} = eva_info{jj,ii}{kk}.Image_VS_Text_precision';
                Text_VS_Image_recall{ii}{jj,kk,:}    = eva_info{jj,ii}{kk}.Text_VS_Image_recall';
                Text_VS_Image_precision{ii}{jj,kk,:} = eva_info{jj,ii}{kk}.Text_VS_Image_precision';

                % Top number Precision
                Image_To_Text_Precision{ii}{jj,kk,:} = eva_info{jj,ii}{kk}.Image_To_Text_Precision;
                Text_To_Image_Precision{ii}{jj,kk,:} = eva_info{jj,ii}{kk}.Text_To_Image_Precision;

                trainT{ii}{jj,kk} = eva_info{jj,ii}{kk}.trainT;
            end
        end
    end

    
    %% Save
    save(result_name,'eva_info','param','loopnbits','hashmethods',...
        'XChunk','XTest','YChunk','YTest','LChunk','LTest',...
        'trainT','Image_VS_Text_MAP','Text_VS_Image_MAP','Image_VS_Text_recall','Image_VS_Text_precision',...
        'Table_ItoT_MAP','Table_TtoI_MAP',...
        'Text_VS_Image_recall','Text_VS_Image_precision','Image_To_Text_Precision','Text_To_Image_Precision','-v7.3');
end