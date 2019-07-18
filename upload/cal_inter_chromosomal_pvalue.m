clear
load('binding2knockout_Reimand2.mat');
load('gene2gene_inter2.mat');
load('genecoordination.mat');
[arow,acol]=size(binding2knockout_refined);
result=zeros(acol,6);
for i=1:acol
    [temprow1,tempcol1]=size(binding2knockout_refined{i}{1});
    [temprow2,tempcol2]=size(binding2knockout_refined{i}{2});
    binding2knockout_refined_ind{i}{1}=zeros(temprow1,1);
    binding2knockout_refined_ind{i}{2}=zeros(temprow2,1);
    binding2knockout_hicpair{i}{1}=[];
    binding2knockout_hicpair{i}{2}=[];
    if temprow1>20 && temprow2>20
        for j=1:temprow1
            ind=0;
            for jj=1:temprow2
                if ind==1
                    break;
                end
                if genecoordination(binding2knockout_refined{i}{1}(j,1),2)~=genecoordination(binding2knockout_refined{i}{2}(jj,1),2) && gene2gene_intra(binding2knockout_refined{i}{1}(j,1),binding2knockout_refined{i}{2}(jj,1))==1
                    ind=1;
                end
            end
            if ind==1
                result(i,1)=result(i,1)+1;
                binding2knockout_refined_ind{i}{1}(j,1)=1;
            end
        end
        for j=1:temprow2
            ind=0;
            for jj=1:temprow1
                if ind==1
                    break;
                end
                if genecoordination(binding2knockout_refined{i}{2}(j,1),2)~=genecoordination(binding2knockout_refined{i}{1}(jj,1),2) && gene2gene_intra(binding2knockout_refined{i}{2}(j,1),binding2knockout_refined{i}{1}(jj,1))==1 
                    ind=1;
                end
            end
            if ind==1
                result(i,2)=result(i,2)+1;
                binding2knockout_refined_ind{i}{2}(j,1)=1;
            end
        end
        result(i,3)=temprow1;
        result(i,4)=temprow2;
        result(i,5)=result(i,1)/result(i,3);
        result(i,6)=result(i,2)/result(i,4);
        ind=0;
        for j=1:temprow1
            for jj=1:temprow2
                if genecoordination(binding2knockout_refined{i}{1}(j,1),2)~=genecoordination(binding2knockout_refined{i}{2}(jj,1),2) && gene2gene_intra(binding2knockout_refined{i}{1}(j,1),binding2knockout_refined{i}{2}(jj,1))==1
                    ind=ind+1;
                    binding2knockout_hicpair{i}{1}(ind,1)=binding2knockout_refined{i}{1}(j,1);
                    binding2knockout_hicpair{i}{1}(ind,2)=binding2knockout_refined{i}{2}(jj,1);
                end
            end
        end
        ind=0;
        for j=1:temprow1
            for jj=1:temprow2
                if genecoordination(binding2knockout_refined{i}{1}(j,1),2)~=genecoordination(binding2knockout_refined{i}{2}(jj,1),2) && gene2gene_intra(binding2knockout_refined{i}{1}(j,1),binding2knockout_refined{i}{2}(jj,1))==1
                    ind=ind+1;
                    binding2knockout_hicpair{i}{2}(ind,1)=binding2knockout_refined{i}{1}(j,1);
                    binding2knockout_hicpair{i}{2}(ind,2)=binding2knockout_refined{i}{2}(jj,1);
                    break;
                end
            end
        end
    end
end
save('binding2knockout_Reimand3_inter.mat','binding2knockout','binding2knockout_refined','result','binding2knockout_refined_ind','binding2knockout_hicpair');
%see varialble 'result' for the 5th column represents the frequency of pairs between TF and its bound genes that show inter-chromosomal colocalization with genes whose expression were affected by knockout of the corresponding TF; the 6th column represents the frequency of TF-KO genes showing inter-chromosomal colocalization with TF bound genes


clear
load('random_binding2knockout_Reimand4.mat');
load('gene2gene_inter2.mat');
load('genecoordination.mat');
[brow,bcol]=size(random_genes);
[arow,acol]=size(random_genes{1});
random_num=bcol;
k=0;
gene2fragments{6767}=[];
for i=1:acol
    result{i}=zeros(random_num,4);
end
for iii=1:bcol
    for i=1:acol
        [temprow1,tempcol1]=size(random_genes{iii}{i}{1});
        [temprow2,tempcol2]=size(random_genes{iii}{i}{2});
        if temprow1>20 && temprow2>20
            for j=1:temprow1
                ind=0;
                for jj=1:temprow2
                    if ind==1
                        break;
                    end
                    if genecoordination(random_genes{iii}{i}{1}(j,1),2)~=genecoordination(random_genes{iii}{i}{2}(jj,1),2) && gene2gene_intra(random_genes{iii}{i}{1}(j,1),random_genes{iii}{i}{2}(jj,1))==1 
                        ind=1;
                    end
                end
                if ind==1
                    result{i}(iii,1)=result{i}(iii,1)+1;
                end
            end
            for j=1:temprow2
                ind=0;
                for jj=1:temprow1
                    if ind==1
                        break;
                    end
                    if genecoordination(random_genes{iii}{i}{1}(jj,1),2)~=genecoordination(random_genes{iii}{i}{2}(j,1),2)  && gene2gene_intra(random_genes{iii}{i}{1}(jj,1),random_genes{iii}{i}{2}(j,1))==1 
                        ind=1;
                    end
                end
                if ind==1
                    result{i}(iii,2)=result{i}(iii,2)+1;
                end
            end
        end
        result{i}(iii,3)=temprow1;
        result{i}(iii,4)=temprow2;
        result{i}(iii,5)=result{i}(iii,1)/result{i}(iii,3);
        result{i}(iii,6)=result{i}(iii,2)/result{i}(iii,4);
    end
end
save('random_binding2knockout_Reimand6_inter.mat','result');
%10,000 randomized experiments


clear
load('binding2knockout_Reimand3_inter.mat');
result1=result;
load('random_binding2knockout_Reimand6_inter.mat');%10,000 randomized experiments
[arow,acol]=size(result1);
result1(arow+1,1)=sum(result1(:,1));
result1(arow+1,2)=sum(result1(:,2));
result1(arow+1,3)=sum(result1(:,3));
result1(arow+1,4)=sum(result1(:,4));
result1(arow+1,5)=result1(arow+1,1)/result1(arow+1,3);
result1(arow+1,6)=result1(arow+1,2)/result1(arow+1,4);
[brow,bcol]=size(result{1});
for i=1:brow
    result{arow+1}(i,1:bcol)=0;
    for j=1:arow
        for jj=1:bcol-2
           result{arow+1}(i,jj)=result{arow+1}(i,jj)+result{j}(i,jj);
        end
    end
    result{arow+1}(i,bcol-1)=result{arow+1}(i,1)/result{arow+1}(i,3);
    result{arow+1}(i,bcol)=result{arow+1}(i,2)/result{arow+1}(i,4);
end
k1=0;
ind=zeros(1,2);
ind1=zeros(arow,2);
for j=1:2
    temp=result{arow+1}(:,j+4);
    temp=sortrows(temp);
    k=brow;
    if result1(arow+1,j+4)<temp(1,1)
        k=0;
    else
        for ii=1:brow
            if result1(arow+1,j+4)>=temp(ii,1)
                k=ii;
            end
        end
    end
    pvalue(1,j)=(brow-k)/brow;
end
%Calculating the P value
%The 1st column represents the P value of pairs between TF and its bound genes that show inter-chromosomal colocalization with genes whose expression were affected by knockout of the corresponding TF; The 2nd column represents the P value of TF-KO genes showing inter-chromosomal colocalization with TF bound genes
