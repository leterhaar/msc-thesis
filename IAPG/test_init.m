addpath('../formulation_SVM');

d = 50;
m = 100;

svm = create_SVM(d,m);

test_sequence = {'IPG', 'IAPG'};