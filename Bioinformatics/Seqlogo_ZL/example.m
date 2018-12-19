%clear
clc

% the program "obtain_lettershapes.m" needs to be run only once
%run('obtain_lettershapes.m')


% this is for generating vector map of sequence logos

% for nuleotide sequence
ntoraa='NT';
colormapping=nuc_colormaping;
weight_thresh=0;


rawseqs = {'CACGTAACATCTC','ACGACGTAACATCTTCT','AAACGTAACATCTCGC'};
alignresults=multialign(rawseqs,'terminalGapAdjust',true);
figure;
maxi=SeqLogo_ZL(alignresults,ntoraa,colormapping,letterclass,weight_thresh);
set(gca,'fontsize',22);
axis([0,size(alignresults,2)+1, 0, maxi]);
xlabel('Position')
ylabel('Bit');

         

% for amino acid sequence
alignedseq_1=repmat(strjoin(Astrs,''),4,1);
alignedseq_2 = {'LSGGQRQRVAIARALAL'; 
              'LSGGEKQRVAIARALMN'; 
              'LSGGQIQRVLLARALAA';
              'LSGGERRRLEIACVLAL'; 
              'FSGGEKKKNELWQMLAL'; 
              'LSGGERRRLEIACVLAL'};
ntoraa='AA';
colormapping=aa_colormaping;
weight_thresh=0;
                   
figure;
subplot(2,1,1)
maxi=SeqLogo_ZL(alignedseq_1,ntoraa,colormapping,letterclass,weight_thresh);
subplot(2,1,2)
maxi=SeqLogo_ZL(alignedseq_2,ntoraa,colormapping,letterclass,weight_thresh);
ylim([0 maxi])

% adjust details of the figure;
