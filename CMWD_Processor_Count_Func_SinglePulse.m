function [data_WJ, data_WOJ] = CMWD_Processor_Count_Func_SinglePulse(Subjects,FilePath)

    data_WJ = [];
    data_WOJ = [];
    
    for FileIndex = 1:length(Subjects)

        fprintf('\n%s\n', Subjects{FileIndex});   

        Data_File = [FilePath 'environment\Rec_NI5105_Updated.bin'];
        fid    = fopen(Data_File,'rb');
        if fid == -1
            fprintf('File cannot open');
        else
          %  fprintf('Environment Data File open\n');
        end

        NumberOfPulses      = 200;
        NoofPlusetoInteg    = 1;  
        RangeSamples        = 1000;
        RangeOffSet         = 17;

        Fs                  = 4000000;                              % Sampling Frequency
        T                   = 1/Fs;                                 % Sample time

        NFFT                = 2^nextpow2(RangeSamples);             % Next power of 2 from length of RangeSamples
        %frequency_axis     = Fs/2*linspace(0,1,NFFT/2+1);
        frequency_axis      = Fs*linspace(0,1,NFFT);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        range               = (frequency_axis*250e-6*3e8)/(2*150e6);                 
        central_frequency   = 66470;                                % centre frequency where gating should be done it depends on distance of target from antenna
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        PulseCounter        = 0;
        RangeSamplesCounter = 1;
        Amplitude           = 9000;

        Environment         = 1;
        Target              = 0;
        WithoutJacketFile   = 0;

        DataI               = zeros(1,RangeSamples); 
        DataQ               = zeros(1,RangeSamples); 
        DataTrigger         = zeros(1,RangeSamples); 

        DataICross          = zeros(NumberOfPulses/2,RangeSamples); 
        DataQCross          = zeros(NumberOfPulses/2,RangeSamples); 

        DataICo             = zeros(NumberOfPulses/2,RangeSamples); 
        DataQCo             = zeros(NumberOfPulses/2,RangeSamples); 

        DataIAndQCross      = zeros(1,NFFT); 
        DataIAndQCo         = zeros(1,NFFT);  

        FFTofDataIAndQCross = zeros(1,NFFT);
        FFTofDataIAndQCo    = zeros(1,NFFT);

        EnvironmentMeanICross = zeros(1,RangeSamples);
        EnvironmentMeanQCross = zeros(1,RangeSamples);

        EnvironmentMeanICo    = zeros(1,RangeSamples);
        EnvironmentMeanQCo    = zeros(1,RangeSamples);

        TriggerPulses       = 0;
        MeanOfTrigger       = -11718; % fix mean of trigger bases on 200 pulses 

        Co_componentM       = [];
        Cross_componentM    = [];

        maxMagOfFFTIAndQCo      =  0;
        maxMagOfFFTIAndQCross   =  0;

        VectorCo     = [];
        VectorCross  = [];

        PeakToAreaRatioCoArray    = [];
        PeakToAreaRatioCrossArray = [];

        MaxIndexCoArray     = [];
        MaxIndexCrossArray  = [];

        while(1)

            DataI               = fread(fid, RangeSamples, 'int16');
            DataQ               = fread(fid, RangeSamples, 'int16');
            DataTrigger         = fread(fid, RangeSamples, 'int16');

            %DataI = filter(ones(1,10)/10,1,DataI);
            %DataQ = filter(ones(1,10)/10,1,DataQ);


            eofstat = feof(fid);
            if eofstat == 1;

                if(Target) % when end of target file reach then Quit the application 
                    fclose(fid); % close the target data file 
                    %fprintf('End of Target data File\n');
                    break;
                end

                %fprintf('End of environment data File\n');
                fclose(fid); % close the environment data file and close of the target data file 

                if(WithoutJacketFile)
                    %WithoutJacketFile = 0;
                    Subject_woj = [Subjects{FileIndex} '_woj'];
                    Data_File = [FilePath Subject_woj '\Rec_NI5105_Updated.bin' ];
                    fid    = fopen(Data_File,'rb');
                    if fid == -1
                     fprintf('File cannot open\n');
                    else
                     %fprintf('Without Jacket Data File open\n');
                    end

                  
    
                else
                    Subject_wj = [Subjects{FileIndex} '_wj'];
                    Data_File = [FilePath Subject_wj '\Rec_NI5105_Updated.bin'];
                    fid    = fopen(Data_File,'rb');
                    if fid == -1
                     fprintf('File cannot open\n');
                    else
                     %fprintf('With Jacket Data File open\n');
                    end
                end

                DataI           = fread(fid, RangeSamples, 'int16');
                DataQ           = fread(fid, RangeSamples, 'int16');
                DataTrigger     = fread(fid, RangeSamples, 'int16');

                Environment         = 0;
                Target              = 1;
                TriggerPulses       = 0;
                PulseCounter        = 0;
                %break;
            end

            %TriggerPulses = [TriggerPulses , DataTrigger'];
            PulseCounter   = PulseCounter + 1;

            if(DataTrigger(1) > MeanOfTrigger)
               Trigger  = 1;
            else
                Trigger = 0;
            end    



            for j = RangeOffSet:RangeSamples-RangeOffSet
                if(Trigger)
                    DataICross(ceil(PulseCounter/2),j) = DataI(j);
                    DataQCross(ceil(PulseCounter/2),j) = DataQ(j);
                else
                    if(PulseCounter == 1 && Trigger == 0)
                        DataICross(ceil(PulseCounter/2),j) = DataI(j);
                        DataQCross(ceil(PulseCounter/2),j) = DataQ(j);
                    else
                        DataICo(ceil(PulseCounter/2),j) = DataI(j);
                        DataQCo(ceil(PulseCounter/2),j) = DataQ(j);
                    end
                end
            end

            if(Environment)

                if(PulseCounter == NumberOfPulses)

                    for rangeindex = RangeOffSet:RangeSamples-RangeOffSet

                        for pulseindex = 1:NumberOfPulses/2

                            DataICrossforMean(pulseindex)  = DataICross(pulseindex,rangeindex);
                            DataQCrossforMean(pulseindex)  = DataQCross(pulseindex,rangeindex);

                            DataICoforMean(pulseindex)  = DataICo(pulseindex,rangeindex);
                            DataQCoforMean(pulseindex)  = DataQCo(pulseindex,rangeindex);

                        end

                            MeanOfPulseICross  = PulseMean(DataICrossforMean,NumberOfPulses/2);
                            EnvironmentMeanICross(rangeindex) = MeanOfPulseICross ;

                            MeanOfPulseQCross  = PulseMean(DataQCrossforMean,NumberOfPulses/2);
                            EnvironmentMeanQCross(rangeindex) = MeanOfPulseQCross;

                            MeanOfPulseICo  = PulseMean(DataICoforMean,NumberOfPulses/2);
                            EnvironmentMeanICo(rangeindex) = MeanOfPulseICo;

                            MeanOfPulseQCo  = PulseMean(DataQCoforMean,NumberOfPulses/2);
                            EnvironmentMeanQCo(rangeindex) = MeanOfPulseQCo;
                    end


                end

            elseif(Target )

                if(PulseCounter == NumberOfPulses)

                    PluseIntegICross = [];
                    PluseIntegQCross = [];
                    PluseIntegICo    = [];
                    PluseIntegQCo    = [];

                    for j = 1:(NumberOfPulses/(2*NoofPlusetoInteg))

                        startIndex = (j-1)*NoofPlusetoInteg + 1;
                        EndIndex   = (j)*NoofPlusetoInteg;

                        PluseIntegICross   = [PluseIntegICross ; (DataICross(startIndex:EndIndex,:))];

                        PluseIntegQCross   = [PluseIntegQCross ; (DataQCross(startIndex:EndIndex,:))];

                        PluseIntegICo      = [PluseIntegICo  ; (DataICo(startIndex:EndIndex,:))];

                        PluseIntegQCo      = [PluseIntegQCo ; (DataQCo(startIndex:EndIndex,:))];
                    end

                    MagOfFFTIAndQCrossArray = [];
                    MagOfFFTIAndQCoArray = [];

                    for pulseindex = 1:(NumberOfPulses/(2*NoofPlusetoInteg))


                        for rangeindex = RangeOffSet:(RangeSamples-RangeOffSet)

                            PluseIntegICross(pulseindex,rangeindex)  = PluseIntegICross(pulseindex,rangeindex)- EnvironmentMeanICross(rangeindex);
                            PluseIntegQCross(pulseindex,rangeindex)  = PluseIntegQCross(pulseindex,rangeindex) - EnvironmentMeanQCross(rangeindex);

                            DataIAndQCross(rangeindex)               = complex(PluseIntegICross(pulseindex,rangeindex) , PluseIntegQCross(pulseindex,rangeindex));

                            PluseIntegICo(pulseindex,rangeindex)     = PluseIntegICo(pulseindex,rangeindex) - EnvironmentMeanICo(rangeindex);
                            PluseIntegQCo(pulseindex,rangeindex)     = PluseIntegQCo(pulseindex,rangeindex) - EnvironmentMeanQCo(rangeindex);

                            DataIAndQCo(rangeindex)                  = complex(PluseIntegICo(pulseindex,rangeindex) , PluseIntegQCo(pulseindex,rangeindex));

                        end

                        FFTofDataIAndQCross    = fft(DataIAndQCross(RangeOffSet+1:RangeSamples-RangeOffSet),NFFT)/(RangeSamples-2*RangeOffSet);
                        FFTofDataIAndQCo       = fft(DataIAndQCo(RangeOffSet+1:RangeSamples-RangeOffSet),NFFT)/(RangeSamples-2*RangeOffSet);

                        MagOfFFTIAndQCross      = abs(FFTofDataIAndQCross);
                        MagOfFFTIAndQCo         = abs(FFTofDataIAndQCo);


                        MagOfFFTIAndQCrossArray = [MagOfFFTIAndQCrossArray;MagOfFFTIAndQCross(15:25)];
                        MagOfFFTIAndQCoArray = [MagOfFFTIAndQCoArray;MagOfFFTIAndQCo(15:25)];

                        index                   = find(frequency_axis >= central_frequency,1);
                        maxMagOfFFTIAndQCo      = 0;
                        maxMagOfFFTIAndQCross   = 0;
                        MaxIndexCo              = 0;
                        MaxIndexCross           = 0;

                        for loopidx = index-2:index+2
                            VectorCo  =  [VectorCo , MagOfFFTIAndQCo(loopidx)]; 
                            if (MagOfFFTIAndQCo(loopidx) > maxMagOfFFTIAndQCo)
                               maxMagOfFFTIAndQCo =  MagOfFFTIAndQCo(loopidx);
                               MaxIndexCo  = loopidx;
                            end
                        end
                        Co_componentM          = [Co_componentM maxMagOfFFTIAndQCo];
                        PeakToAreaRatioCo      =  maxMagOfFFTIAndQCo/mean(VectorCo);
                        PeakToAreaRatioCoArray = [PeakToAreaRatioCoArray PeakToAreaRatioCo];
                        MaxIndexCoArray        = [MaxIndexCoArray MaxIndexCo];    


                        for loopidx = index-2:index+2
                            VectorCross  =  [VectorCross, MagOfFFTIAndQCross(loopidx)]; 
                            if (MagOfFFTIAndQCross(loopidx) > maxMagOfFFTIAndQCross)
                               maxMagOfFFTIAndQCross =  MagOfFFTIAndQCross(loopidx);
                               MaxIndexCross  = loopidx;
                            end
                        end
                        Cross_componentM          = [Cross_componentM maxMagOfFFTIAndQCross];
                        PeakToAreaRatioCross      =  maxMagOfFFTIAndQCross/mean(VectorCross);
                        PeakToAreaRatioCrossArray = [PeakToAreaRatioCrossArray PeakToAreaRatioCross];
                        MaxIndexCrossArray        = [MaxIndexCrossArray MaxIndexCross];

                        VectorCo    = [];
                        VectorCross = [];



                    end  % end of for pulseindex = 1:NumberOfPulses/2

                    stdCo                         = std(Co_componentM);
                    stdCross                      = std(Cross_componentM);

                    Co_componentMArray            = Co_componentM;
                    Cross_componentMArray            = Cross_componentM;

                    Co_componentM                 = mean(Co_componentM);
                    Cross_componentM              = mean(Cross_componentM);

                    MaxIndexCoArrayMedian         = mode(MaxIndexCoArray);
                    MaxIndexCrossArrayMedian      = mode(MaxIndexCrossArray);




                    fprintf('Co Mean = %f, Co Std = %f , Cross  Mean = %f, Cross Std=%f \n',Co_componentM,stdCo,Cross_componentM,stdCross); 


                    if(WithoutJacketFile == 1)


                            dataForFFT_Cross_WOJ = [MagOfFFTIAndQCrossArray];
                            Cross_Count_Bw_50_65 = length(find(dataForFFT_Cross_WOJ(:,4)>= 48));
                            Cross_Count_Bw_50_65 = Cross_Count_Bw_50_65+length(find(dataForFFT_Cross_WOJ(:,5)>= 48));
                            feature_vec = [Cross_Count_Bw_50_65];
                            data_WOJ = [data_WOJ; feature_vec];
                    else

                        
                            dataForFFT_Cross_WJ = [MagOfFFTIAndQCrossArray];
                            Cross_Count_Bw_50_65 = length(find(dataForFFT_Cross_WJ(:,4)>= 48));
                            Cross_Count_Bw_50_65 = Cross_Count_Bw_50_65+length(find(dataForFFT_Cross_WJ(:,5)>= 48));
                            feature_vec = [Cross_Count_Bw_50_65];
                            data_WJ = [data_WJ; feature_vec];
                        
                    end


                    if(WithoutJacketFile  == 1)

                        break;

                    end

                    Target              = 0; % make target is equal to zero, for the opeing of faheem without jacket file 
                    WithoutJacketFile   = 1; % make WithoutJacketFile is equal to 1, for the opeing of faheem without jacket file 
                    stdCo               = 0;
                    stdCross            = 0;
                    Co_componentM       = [];
                    Cross_componentM    = [];
                   PeakToAreaRatioCoArray = [];
                   PeakToAreaRatioCrossArray = [];
                end  % end of if(PulseCounter == NumberOfPulses)

            else
                %fprintf('Else of Environment&Target \n');
            end

            %fprintf('testing \n');
        end
        %fclose(MLSubjectFile);
    end

end