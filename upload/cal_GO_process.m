clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('GO_process.mat');
[arow,acol]=size(result);
[brow,bcol]=size(gene2gene{1});
num5=0;
num6=0;
num7=0;
num8=0;
result2=zeros(arow,4);
for i=1:arow
    k1=0;
    k2=0;
    temp_ind1=[];
    temp_ind2=[];
    [temprow,tempcol]=size(binding2knockout_refined_ind{i}{2});
    for j=1:temprow
        if binding2knockout_refined_ind{i}{2}(j,1)==1 || binding2knockout_refined_ind1{i}{2}(j,1)==1
            k1=k1+1;
            temp_ind1(k1,1)=binding2knockout_refined{i}{2}(j,1);
            geneid{i}(k1,1)=binding2knockout_refined{i}{2}(j,1);
        else
            k2=k2+1;
            temp_ind2(k2,1)=binding2knockout_refined{i}{2}(j,1);
        end
    end
    [crow,ccol]=size(binding2knockout{i}{1});
     for ii=1:k1
        num6=num6+1;
        result2(i,2)=result2(i,2)+1;
        for jj=1:crow
            if gene2gene{1}(temp_ind1(ii,1),binding2knockout{i}{1}(jj,1))~=0
                num5=num5+1;
                result2(i,1)=result2(i,1)+1;
                break;
            end
        end
    end
    for ii=1:k2
        num8=num8+1;
        result2(i,4)=result2(i,4)+1;
        for jj=1:crow
            if gene2gene{1}(temp_ind2(ii,1),binding2knockout{i}{1}(jj,1))~=0
                num7=num7+1;
                result2(i,3)=result2(i,3)+1;
                break;
            end
        end
    end
end
result1(1,1)=num5/num6;
result1(1,2)=num7/num8;
%Values in 'result1' correspond to the frequencies of all pairs between TFs and their SAR genes (or TF non-SAR genes) that are involved in the same biological process with genes bound by the corresponding TFs
