clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('TF_binding_genes_ypd.mat');
load('GO_component.mat');
a=load('test.txt');
b=load('test1.txt');
[arow,acol]=size(TF_binding_genes);
[crow,ccol]=size(binding2knockout_refined);
TF_binding_genes1=zeros(6767,acol);
for i=1:acol
    [brow,bcol]=size(TF_binding_genes{i});
    for j=1:brow
        TF_binding_genes1(TF_binding_genes{i}(j,1),i)=1;
    end
end
k1=0;
k2=0;
for i=1:ccol
    [brow,bcol]=size(binding2knockout_refined{i}{2});
    k3=0;
    tempid=[];
    for ii=1:brow
        if binding2knockout_refined_ind{i}{2}(ii,1)==1 || binding2knockout_refined_ind1{i}{2}(ii,1)==1
            k3=k3+1;
            tempid(k3,1)=binding2knockout_refined{i}{2}(ii,1);
        end
    end
    for jj=1:k3
        for j=2:acol
            if TF_binding_genes1(tempid(jj,1),j)==1 && gene2gene{1}(a(j,1),b(i,1))~=0
                k1=k1+1;
                break;
            end
        end
    end
    k2=k2+k3;
end
ratio=k1/k2;
%The resulting ratio represents frequency of pairs between TFs and their SAR genes showing cellular co-component between their SAR TFs and binding TFs of TF SAR genes



clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('TF_binding_genes_ypd.mat');
load('GO_component.mat');
a=load('test.txt');
b=load('test1.txt');
[arow,acol]=size(TF_binding_genes);
[crow,ccol]=size(binding2knockout_refined);
TF_binding_genes1=zeros(6767,acol);
for i=1:acol
    [brow,bcol]=size(TF_binding_genes{i});
    for j=1:brow
        TF_binding_genes1(TF_binding_genes{i}(j,1),i)=1;
    end
end
random_num=10000;
for iii=1:random_num
    k1=0;
    k2=0;
    for i=1:ccol
        [brow,bcol]=size(binding2knockout_refined{i}{2});
        k3=0;
        tempid=[];
        for ii=1:brow
            if binding2knockout_refined_ind{i}{2}(ii,1)==1 || binding2knockout_refined_ind1{i}{2}(ii,1)==1
                k3=k3+1;
                tempid(k3,1)=ceil(6767*rand);
            end
        end
        for jj=1:k3
            for j=2:acol
                if TF_binding_genes1(tempid(jj,1),j)==1 && gene2gene{1}(a(j,1),b(i,1))~=0
                    k1=k1+1;
                    break;
                end
            end
        end
        k2=k2+k3;
    end
    ratio(iii,1)=k1/k2;
end
save('random_componet1.mat','ratio');
%10,000 randomized experiments