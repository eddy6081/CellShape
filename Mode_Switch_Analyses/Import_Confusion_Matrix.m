function [CM] = Import_Confusion_Matrix()
    %{
    Import confusion matrix saved from All_SVM.py
    %}
    %be all we need.
    disp('Please select csv file with confusion matrix saved from All_SVM.py:   ')
    [fil,pth]=uigetfile('*.csv');
    
    T = readtable(fullfile(pth,fil));
    CM = T{:,:};
end