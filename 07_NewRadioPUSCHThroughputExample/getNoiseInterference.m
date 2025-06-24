function [NI] = getNoiseInterference(decbits, encodeULSCH, harqEntity, pusch, carrier, NTxAnts, rxGrid, H_DMRS)
    
    fprintf("harqEntity.RedundancyVersion: %d\n", harqEntity.RedundancyVersion);

    [puschIndices, puschIndicesInfo] = nrPUSCHIndices(carrier, pusch);

    setTransportBlock(encodeULSCH, decbits, harqEntity.HARQProcessID)
    codedTrBlock = encodeULSCH(pusch.Modulation, pusch.NumLayers, puschIndicesInfo.G, harqEntity.RedundancyVersion);
    puschGrid = nrResourceGrid(carrier, NTxAnts);
    puschSymbols = nrPUSCH(carrier,pusch,codedTrBlock);

    if (strcmpi(pusch.TransmissionScheme,'codebook'))
        F = eye(pusch.NumAntennaPorts, NTxAnts);
    else
        F = eye(pusch.NumLayers, NTxAnts);
    end

    [~,puschAntIndices] = nrExtractResources(puschIndices, puschGrid);
    puschGrid(puschAntIndices) = puschSymbols * F;

    dmrsSymbols = nrPUSCHDMRS(carrier, pusch);
    dmrsIndices = nrPUSCHDMRSIndices(carrier,pusch);

    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),puschGrid);
        puschGrid(dmrsAntIndices) = puschGrid(dmrsAntIndices) + dmrsSymbols(:,p) * F(p,:);
    end

    NI = rxGrid - puschGrid.*H_DMRS;
end
