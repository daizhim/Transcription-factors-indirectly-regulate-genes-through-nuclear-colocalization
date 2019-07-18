clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('TF_noboundgene_MacIsaac.mat');
clear result
[arow,acol]=size(binding2knockout_refined);
[drow,dcol]=size(geneid2);
k1=0;
k11=0;
k111=0;
ind=zeros(acol,3);
for i=1:acol
    [brow,bcol]=size(geneid1{i});
    if brow~=0
        [crow,ccol]=size(binding2knockout_refined{i}{1});
        k=0;
        for ii=1:crow
            if geneid2(binding2knockout_refined{i}{1}(ii,1),i)==1
                k=k+1;
            end
        end
        [crow,ccol]=size(binding2knockout_refined{i}{2});
        k=0;
        k2=0;
        k4=0;
        k3=0;
        for ii=1:crow
            if binding2knockout_refined_ind{i}{2}(ii,1)==1 || binding2knockout_refined_ind{i}{2}(ii,1)==1
                k2=k2+1;
                if geneid2(binding2knockout_refined{i}{2}(ii,1),i)==1
                    k=k+1;
                end
            else
                k4=k4+1;
                if geneid2(binding2knockout_refined{i}{2}(ii,1),i)==1
                    k3=k3+1;
                end
            end
        end
        result(i,1)=k;
        result(i,2)=k2-k;
        result(i,3)=k3;
        result(i,4)=k4-k3;
    end
end
result1=sum(result);
ave(1,1)=result1(1,1)/(result1(1,1)+result1(1,2));
ave(1,2)=result1(1,3)/(result1(1,4)+result1(1,3));
%The ave are values that correspond to the frequencies of pairs between all TFs and their SAR genes (or TF non-SAR genes) that have DNA motifs of the corresponding TFs in promoter regions.
