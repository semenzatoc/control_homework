function [patient] = patient_parameters(subj)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
Patients=[1   74   164   88   1  32.7  60   98  2   5.2     85;
          2   67   161   69   1  26.6  53   96  1   5.8     81;
          3   75   176   101  1  32.6  69   100 4   6.1     79;
          4   69   173   97   1  32.4  67   91  3   5.6     86;
          5   45   171   64   1  21.9  52   94  2   6       82;
          6   57   182   80   1  24.2  62   89  1   6.2     83;
          7   74   155   55   1  22.9  44   95  2   5.1     89;
          8   71   172   78   1  26.4  60   86  1   5.9     80;
          9   65   176   77   1  24.9  60   100 6   4.8     84;
          10  72   192   73   1  19.8  62   98  3   6.3     80;
          11  69   168   84   1  29.8  60   90  2   5.5     87;
          12  60   190   92   1  25.5  71   94  3   5.6     79;
          13  61   177   81   1  25.9  62   96  4   6.0     82;
          14  54   173   86   1  27.5  63   93  4   5.1     83;
          15  71   172   83   1  28.1  62   87  2   4.9     78;
          16  53   186   114  1  33    77   97  3   4.5     81;
          17  72   162   87   1  33.2  59   85  1   6.3     77;
          18  61   182   93   1  28.1  69   100 5   5.8     80;
          19  70   167   77   1  27.6  58   97  2   5.3     82;
          20  69   168   82   1  29.1  60   92  1   5.0     84;
          21  69   158   81   1  32.4  55   98  4   5.6     88;
          22  60   165   85   1  31.2  60   99  4   6.3     86;
          23  70   173   69   1  23.1  56   88  1   6.1     80;
          24  56   186   99   1  28.6  73   95  2   4.9     80];     
[m, n]=size(Patients); 
if subj > m
    error('Subject must be maximum 24')
end
if subj < 0
    error('Subject number must be a positive integer or zero (for the average subject)')
end
if subj>0
    patient.age = Patients(subj,2); 
    patient.height = Patients(subj,3); 
    patient.weight = Patients(subj,4); 
    patient.Gender = Patients(subj,5); 
    patient.BMI = Patients(subj,6); 
    patient.lbm = Patients(subj,7);
    patient.E0 = Patients(subj,8);
    patient.Emax= patient.E0 - Patients(subj,9);
%     patient.Cobasis= Patients(subj,10); % 5L/min
%     patient.MAPbasis= Patients(subj,11); % 80mmHg 
else 
    
    patient.age = mean(Patients(:,2)); 
    patient.height = mean(Patients(:,3)); 
    patient.weight = mean(Patients(:,4)); 
    patient.Gender = 1; 
    patient.BMI = mean(Patients(:,6)); 
    patient.lbm = mean(Patients(:,7));
    patient.E0 = mean(Patients(:,8));
    patient.Emax= mean(patient.E0 - Patients(:,9));
%     patient.Cobasis= mean(Patients(:,10)); % 5L/min
%     patient.MAPbasis= mean(Patients(:,11)); % 80mmHg 
end
end

