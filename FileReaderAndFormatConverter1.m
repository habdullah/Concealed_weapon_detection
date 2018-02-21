clear all;
close all;
clc;
wj = 0;


%Subjects = {'hassan','adeel1','Ameen1','amjad1','Arqam1','Asif1','Azaman1','bhatti1','faheem1','furkan1','george1','gohas1','grizwan1','hassan1','kashif1','naseem1','naveed1','peter1','rafiq1','rizwan1','sajjad1','sakim1','tariq1','umer1','usman1','waqas1','waqash1'};
Subjects ={'abid','abid1','abid2','abidz','abidz1','abidz2','abida','abida1','abida2','adil','adil1','adil2','ahsan','ahsan1','ahsan2','ali','ali1','ali2','alishan','alishan1','alishan2','aliwaqas','aliwaqas1','aliwaqas2','ameen','ameen1','ameen2','amir','amir1','amir2','amjad','amjad1','amjad2','aqayyum','aqayyum1','aqayyum2','asif','asif1','asif2','azhar','azhar1','azhar2','furqan','furqan1','furqan2','gm','gm1','gm2','hassan','hassan1','hassan2','jahanzaib','jahanzaib1','jahanzaib2','khalid','khalid1','khalid2','khan','khan1','khan2', 'khurshid','khurshid1','khurshid2','kirdad','kirdad1','kirdad2','majid','majid1','majid2','mushtaq','mushtaq1','mushtaq2','naseem','naseem1','naseem2','nawaz','nawaz1','nawaz2','noman','noman1','noman2','sakim','sakim1','sakim2','shahid','shahid1','shahid2','shakir','shakir1','shakir2','tariq','tariq1','tariq2','usman','usman1','usman2','waqasI','waqasI1','waqasI2','waseem','waseem1','waseem2','zahid','zahid1','zahid2','zardad','zardad1','zardad2'};

for i = 1:length(Subjects)
    

       if(wj == 1)
          Subject_File = [Subjects{i} '_wj' ];
       else
           Subject_File = [Subjects{i} '_woj' ];
       end



    Data_FileI = ['C:\Data\17-11-2016\' Subject_File '\Rec_NI5105_ch1.bin'];
    fidI    = fopen(Data_FileI,'rb');
    if fidI == -1
     fprintf('File cannot open\n');
    else
     fprintf('File Rec_NI5105_ch1.bin open\n');
    end

    Data_FileQ = ['C:\Data\17-11-2016\' Subject_File '\Rec_NI5105_ch2.bin'];
    fidQ    = fopen(Data_FileQ,'rb');
    if fidQ == -1
     fprintf('File cannot open\n');
     else
     fprintf('File Rec_NI5105_ch2.bin open\n');
    end

    Data_FileT = ['C:\Data\17-11-2016\' Subject_File '\Rec_NI5105_ch3.bin'];
    fidT    = fopen(Data_FileT,'rb');
    if fidT == -1
     fprintf('File cannot open\n');
     else
     fprintf('File Rec_NI5105_ch3.bin open\n');
    end

    Data_File = ['C:\Data\17-11-2016\' Subject_File '\Rec_NI5105_Updated.bin'];
    fid    = fopen(Data_File,'wb');
    if fid == -1
        fprintf('File cannot open');
    end


    NumberOfPulses      = 200;
    RangeSamples        = 1000;
    PulseCounter        = 1;
    RangeSamplesCounter = 1;


    DataI               = zeros(NumberOfPulses,RangeSamples); 
    DataQ               = zeros(NumberOfPulses,RangeSamples); 
    DataT               = zeros(NumberOfPulses,RangeSamples); 

    DataReadI = fread(fidI,inf,'int16');
    DataReadQ = fread(fidQ,inf,'int16');
    DataReadT = fread(fidT,inf,'int16');


    for i = 1:length(DataReadI)

        if ( mod(i,RangeSamples) ~= 0)

            DataI(PulseCounter,RangeSamplesCounter) = DataReadI(i);
            DataQ(PulseCounter,RangeSamplesCounter) = DataReadQ(i);
            DataT(PulseCounter,RangeSamplesCounter) = DataReadT(i);

            RangeSamplesCounter = RangeSamplesCounter+ 1;
        else

            DataI(PulseCounter,RangeSamplesCounter) = DataReadI(i);
            DataQ(PulseCounter,RangeSamplesCounter) = DataReadQ(i);
            DataT(PulseCounter,RangeSamplesCounter) = DataReadT(i);

    %         subplot(311);plot(DataI(PulseCounter,:));
    %         subplot(312);plot(DataQ(PulseCounter,:));
    %         subplot(313);plot(DataT(PulseCounter,:));

            fwrite(fid, DataI(PulseCounter,:), 'int16');
            fwrite(fid, DataQ(PulseCounter,:), 'int16');
            fwrite(fid, DataT(PulseCounter,:), 'int16');

            PulseCounter        = PulseCounter + 1;
            RangeSamplesCounter = 1;

    %         pause(1);
        end
    end


        fclose(fidI);
        fclose(fidQ);
        fclose(fidT);
        fclose(fid);

end

fprintf('Data Format are convert\n');
fprintf('End of Application\n');

