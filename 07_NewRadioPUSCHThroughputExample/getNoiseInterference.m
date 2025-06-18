function [H] = getNoiseInterference(decbits, encodeULSCH, harqProcessId, pusch, carrier, NTxAnts, rxGrid)
    
    [puschIndices, puschIndicesInfo] = nrPUSCHIndices(carrier, pusch);

    setTransportBlock(encodeULSCH, decbits, harqProcessId)
    codedTrBlock = encodeULSCH(pusch.Modulation, pusch.NumLayers, puschIndicesInfo.G, harqProcessId);
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

    H = puschGrid./rxGrid;
end
