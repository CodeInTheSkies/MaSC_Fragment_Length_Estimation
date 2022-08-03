% -------------------------------------------------------------------------------
% This is the MATLAB prototype written from conception of the idea. 
% This code is a bit raw, which was then used to refine and 
% develop, test, and finalize the Perl implementation
% Therefore, the inline comments you will see are a bit casual and unrefined.
% -------------------------------------------------------------------------------
% This MATLAB code performs the correlation over the whole genome by appending the
% chromosomes one after the other. The read starts and the read ends are loaded 
% from MAT files for the + and - strands, respectively, and then a simple procedure 
% using the definition of correlation is implemented
% ****
% The code also loads the mappability map for all the chromosomes and then accounts
% for mappability

function [CurShft,RdLength] = MaSC_RAW_Implementation(Dflt_MAT_Fldr, MappabFldr, FigsSvDir, PNGSvDir, fnam)

datestr(clock)

disp(' ')
disp('===================================')
disp(['FULL GENOME - With Mappability' fnam])
disp('===================================')
disp(' ')

% Initialization
D = -50:450;
CorrSig = [];
ChrGaps = 2000;
RdsRejectVec = [];
ToTNoOrgRds = 0;

% Getting the list of chromosomes
filLst = dir([Dflt_MAT_Fldr '/*.mat']);
chrList = {};
for jis1 = 1:length(filLst)
    if ~isempty(strfind(filLst(jis1).name, '_Plus'))
        FLNm = filLst(jis1).name;
        part = FLNm(1:strfind(FLNm, '_Plus')-1);
        ind = strfind(part,'_');
        chrnam = part(ind(end)+1:end);
        chrList = [chrList {chrnam}];
    end
    if ~isempty(strfind(filLst(jis1).name, '_Minus'))
        FLNm = filLst(jis1).name;
        part = FLNm(1:strfind(FLNm, '_Minus')-1);
        ind = strfind(part,'_');
        chrnam = part(ind(end)+1:end);
        chrList = [chrList {chrnam}];
    end
end
chrList = unique(chrList);

disp('Getting Mappability data ...')
% Positive strand map data
for jisw = 1:length(chrList)
    CurchrNam = chrList{jisw};
    MapDatt = Load_Mappability_UsingMAT(MappabFldr, CurchrNam);
    Mapab.(chrList{jisw}) = MapDatt;
end

% NOTE: Here, although we have both the positive and the negative strand 25mer mappability
% data available as separate files obtained from the UniqueOme website, we are using only the
% positive data, as the code we wrote earlier is such that we need only the positive
% data; On reviewing the earlier code, we understand that this is because we consider
% only the read starts and not the read ends. Hence, the map is the same for both the
% strands. In case you need to use the negative map, then uncomment the next few lines
% to read the negative map. 
% You should also then uncomment the line "Mapab.(CurchrNam) = MapDatt"
% above so that the positive and the negative map data are properly labeled in the
% structure.

% Retrieving the original reads from the chromosomes and storing in the workspace
disp('Retrieving the original reads')
for jisw = 1:length(chrList)
    CurchrNam = chrList{jisw};
    % Call the function to obtain the sorted Plus and Minus reads for the chromosome
    [CurSrtdPlusRdStrts, CurSrtdPlusRdEnds, CurSrtdMinusRdStrts, CurSrtdMinusRdEnds] = FullChrShftEst_GettheReads_UsingMAT(Dflt_MAT_Fldr, CurchrNam);
    OrigReads.(CurchrNam).PlusRdStrts = CurSrtdPlusRdStrts;
    OrigReads.(CurchrNam).PlusRdEnds = CurSrtdPlusRdEnds;
    OrigReads.(CurchrNam).MinusRdStrts = CurSrtdMinusRdStrts;
    OrigReads.(CurchrNam).MinusRdEnds = CurSrtdMinusRdEnds;
    ToTNoOrgRds = ToTNoOrgRds + length(CurSrtdPlusRdStrts) + length(CurSrtdMinusRdStrts);
end
% Computing the Read length - Used for adjusting the final shift estimate
RdLength = OrigReads.chr1.PlusRdEnds(1) - OrigReads.chr1.PlusRdStrts(1) + 1;
disp(['Calculated Read length = ' num2str(RdLength)])

%%% NOTE: Please note that for chrY, the mappability starts at 0, which gives an
%%% index error in Matlab!  So, just for convenience, we have manually set that 0 to
%%% 1. Of course, we lose that one mappable base pair, but this should not affect the
%%% results in any noticeable manner.

% Now comes the actual computation of the corr. coeff.!

disp('Computing Corr. Coeff...')
cntt = 1;
for d = D
    PlusVec = [];
    MinusVec = [];
    GenSz = 0;
    disp(' '), disp(['Shift ' num2str(d)])
    
    %     Getting the new set of + and - Reads as well as the new Genome size for this
    %     particular shift
    for jisw = 1:length(chrList)
        CurchrNam = chrList{jisw};
        fprintf('%s', [CurchrNam ' '])
        TempRngStrts = Mapab.(CurchrNam).RngStrts;
        TempRngEnds  = Mapab.(CurchrNam).RngEnds;
        GenSz = GenSz + Intersect_Ranges_with_Ranges(TempRngStrts, TempRngEnds, TempRngStrts-d, TempRngEnds-d);   

        % Filtering the original reads for its own strand mappability (we
        % have to do this because of the difference in the UCSC mappability
        % and the mappability we have for our reads. This has to do with
        % the differences in the ELAND criteria and the UCSC criteria
        % Filtering the original set of reads against the Mappability data 

        % NOTE: From here, we are taking only the RdStrts for both + and - strands
        % So, we have to adjust the estimated shift values to accommodate the read
        % length i.e., add RdLength-1 to the shift estimate. Instead, if you take the
        % ends of the - reads, then you needn't do this adjustment. However, you would
        % have to do other adjustments like adjusting for the fact that the read ends
        % are one-based. That would complicate the implementation

        if cntt == 1
            FiltReads.(CurchrNam).PlusRdStrts = ...
                Intersect_NumBrs_with_Ranges(OrigReads.(CurchrNam).PlusRdStrts, TempRngStrts, TempRngEnds);
            FiltReads.(CurchrNam).MinusRdStrts = ...
                Intersect_NumBrs_with_Ranges(OrigReads.(CurchrNam).MinusRdStrts, TempRngStrts, TempRngEnds);
        end
        CurSrtdPlusRdStrts  = ... % New + reads
            Intersect_NumBrs_with_Ranges(FiltReads.(CurchrNam).PlusRdStrts, TempRngStrts-d, TempRngEnds-d);        
        CurSrtdMinusRdStrts = ... % New - reads 
            d + Intersect_NumBrs_with_Ranges(FiltReads.(CurchrNam).MinusRdStrts-d, TempRngStrts, TempRngEnds); 
% (NOTE: in the above line, 'd' must be added to the result to compensate. i.e., to obtain the original reads after filtering)
        % Append the read starts to the long vector
        if jisw == 1
            PlusVec = CurSrtdPlusRdStrts;
%             MinusVec = CurSrtdMinusRdEnds;
            MinusVec = CurSrtdMinusRdStrts;
        else
            ep = PlusVec(end);
            en = MinusVec(end);
            if ep > en
                cn = ep-en;
                cp = 0;
            elseif en > ep
                cp = en-ep;
                cn = 0;
            else
                cp = 0;
                cn = 0;
            end
            PlusVec = [PlusVec ep+cp+ChrGaps+CurSrtdPlusRdStrts];       % With      correction factor
            MinusVec = [MinusVec en+cn+ChrGaps+CurSrtdMinusRdStrts];    % With      correction factor            
        end
    end
    
    RdsRejectVec = [RdsRejectVec ToTNoOrgRds-length(PlusVec)-length(MinusVec)];
    
    %     Computing Means and Std deviations
    PMean = length(PlusVec)/GenSz;
    MMean = length(MinusVec)/GenSz;
    Pstd = sqrt(PMean - PMean^2);
    Mstd = sqrt(MMean - MMean^2);
    
    %     Computing Corr. coeff.
    CorrSig(end+1) = (length(intersect(PlusVec, MinusVec-d))/GenSz - PMean*MMean)/(Pstd*Mstd);
    cntt = cntt+1;
end

% Smoothing the correlation signal CorrSig using a simple moving avg. filter
% Wrong Smoothing algorithm - Signal was offset to the left
% MVaN = 10;  % No. of points for the moving average
% CorrSig = filter((1/MVaN)*ones(1,MVaN), 1, CorrSig(:)');
% Correct Smoothing algorithm
for jis = 1:length(CorrSig)
    SmCorrSig(jis) = mean(CorrSig(max(1,jis-15):min(jis+15,length(CorrSig))));
end

[MaxVlu, MaxIdX] = max(SmCorrSig);
CurShft = D(MaxIdX);

%  Correcting for the read length
CurShft = CurShft + RdLength - 1;

figure,
plot(D, SmCorrSig, 'b', 'linewidth', 2), hold on
plot(D, CorrSig, 'r', 'linewidth', .5),
title({['WithMAP CorrPlot-' strrep(fnam,'_','-') '- Full Genome']; ['Final Estd FL: ' num2str(CurShft)]})

saveas(gcf, [FigsSvDir '/WithMAP_CorrPlot_' strrep(fnam,'-','_') '_FullGenome.fig'], 'fig'),
% save([FigsSvDir '/CorrPlot_' strrep(fnam,'-','_') '_FullGenome.mat'], 'SmCorrSig', 'CorrSig');

% saveas(gcf, [PNGSvDir '/WithMAP_CorrPlot_' strrep(fnam,'-','_') '_FullGenome.png'], 'png'),

disp(['Full Genome / Corrected Shift = ' num2str(CurShft)])
disp('*************************')
disp(' ')

% figure,
% plot(D, RdsRejectVec, '*')
% title(['No. of Rds Rejected vs Shift for ' strrep(fnam,'-','_') ' | Orig No. Rds ' num2str(ToTNoOrgRds)])
% saveas(gcf, [FigsSvDir '/Rds_Rejected_' strrep(fnam,'-','_') '_FullGenome.fig'], 'fig'),

datestr(clock)

close all






