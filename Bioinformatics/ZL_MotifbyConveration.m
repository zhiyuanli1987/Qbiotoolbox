function [condiH_eachsite,totalH,meandist_eachsite]=ZL_MotifbyConveration(alignedseq,ntoraa,distmofelements,printornot)
% this program is for identify nt or aa based distance between pre-aligned
% sequences, and find out motifs if necessary

% start with aligned sequences, check if they are aligned well
seqnum=size(alignedseq,1);
seqlen=size(alignedseq,2);
if seqnum>1 && seqlen>1
    if printornot
        fprintf('seq num (%d), seq len (%d)\n',seqnum,seqlen);
    end
    
    % aa or nt Distance on each site:
    % translate the sequence into numbers
    if strcmp(ntoraa,'nt')
        alignedint=nt2int(alignedseq);
        voidnum=nt2int('-');
    else
        alignedint=aa2int(alignedseq);
        voidnum=aa2int('-');
    end
    
    %distance scoring matrix
    if isempty(distmofelements) % 1 for not the same, 0 for the same
        distmofelements=ones(max(alignedint(:)));
        for i=1:size(distmofelements,1)
            distmofelements(i,i)=0;
        end
        distmofelements(voidnum,voidnum)=1;
    end
    
    %mean distance along each site & mutual information along each site
    frequency_aa=tabulate(alignedint(:));
    frequency_aa=frequency_aa(:,3)/100;
    frequency_aa(frequency_aa==0)=[];
    totalH=-1*sum(frequency_aa.*log2(frequency_aa),1);
    
    meandist_eachsite=nan*zeros(1,seqlen);
    condiH_eachsite=nan*zeros(1,seqlen);
    for s=1:seqlen
        if sum(alignedint(:,s)~=voidnum)>0 %  not all are void
        tabfreq=tabulate(alignedint(:,s));
        tabfreq(tabfreq(:,2)==0,:)=[];
        freq=tabfreq(:,3)/100;
        condiH_eachsite(s)=-1*sum(freq.*log2(freq),1);
        
        matrix_inorder_pairfreq=zeros(max(tabfreq(:,1)));        
        matrix_inorder_pairfreq(tabfreq(:,1),tabfreq(:,1))=freq*freq';
        meandist_eachsite(s)=sum(sum(distmofelements(1:size(matrix_inorder_pairfreq,1),1:size(matrix_inorder_pairfreq,1)).*matrix_inorder_pairfreq));
        end
    end
    
    %mutual information along each site
else
    perror('wrong number of sequences')
end
