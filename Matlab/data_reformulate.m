function InputFile=data_reformulate(data,label)

   data_dumple=[label,data];

    fcost=fopen('data\dumple.txt', 'w'); 
    fclose(fcost);
    fcost=fopen('data\dumple.txt', 'a'); 
    for iii=1:size(data_dumple,1)
        for jjj=1:size(data_dumple,2)
             fprintf(fcost,'%5.12f ',data_dumple(iii,jjj));
        end
         fprintf(fcost,'\n');
    end
    fclose(fcost);
    InputFile = 'data\dumple.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%end of dumple
% 