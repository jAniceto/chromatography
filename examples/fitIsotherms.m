function fitIsotherms()

Cq = [  0	0
        0.142	0.459
        0.153	0.374
        0.649	0.599
        0.672	0.572
        1.400	0.737
        2.073	0.820
        2.162	0.796
        2.945	0.809
        3.530	0.775
        3.596	0.827
        3.696	0.779
        3.729	0.799];
   
C = Cq(:,1);
q = Cq(:,2);

addpath('../');  % add parent folder to path
fitIsotherm(C, q, [1 2], 'langmuir', 'ls')