function maxi=SeqLogo_ZL(alignedseq,ntoraa,colormapping,letterclass,weight_thresh)

% input:
% alignedseq: nucleotide or amino acid sequences that are already aligned, should be a group of sequences in same length
% ntoraa:'NT'(for nucleotide sequence) or 'AA' (for amino acid sequence)
% colormapping: a matlab map from letter to color, can be customized
% letteclass: a structure containing shapes description for letters, need
% to run the program "" to generate
% weight_thresh: if the weight of a letter is smaller than the threshold,
% it will not be shown in the plot

% return:
% a vector plot for manipulation
% maximal height for letters in the plot (for manipulating the plot)

hold on;

% this progra ultilized the built-in program seqlogo in the bioinformatics
% toolbox of Matlab
weightlogo= seqlogo(alignedseq,'Alphabet',ntoraa,'Displaylogo','false');

logoletters=weightlogo{1};
logoweights=weightlogo{2};

maxi=0;
for s=1:size(logoweights,2)
    logoweight_thiss=logoweights(:,s);
    
    [value,order]=sort(logoweight_thiss);
    letterstart=0;
    for i=1:length(order)
        
        lettertodisp=logoletters(order(i));
        letterheight=value(i);
        if letterheight>weight_thresh
            % plot the letter
            letterid=letterclass.map(lettertodisp);
            letterstr=letterclass.strs{letterid};
            letterxcenter=letterclass.centerx(letterid);
            lettercolor=colormapping(letterstr);
            curvenum=letterclass.curvenum(letterid);
            letter_xy=letterclass.curves_xy{letterid};
            for c=1:curvenum
                if c==1
                    tofillcolor=lettercolor;
                else
                    tofillcolor=[1,1,1];
                end
                xs=letter_xy{c}(:,1);
                ys=letter_xy{c}(:,2);
                % adjust x and y
                xs=xs-letterxcenter+s;
                ys=letterstart+letterheight*ys;
                fh=fill(xs,ys,tofillcolor);
                fh.LineStyle='none';
            end
            
            
            letterstart=letterstart+letterheight;
        end
        
    end
    maxi=max(maxi,letterstart);
    
end