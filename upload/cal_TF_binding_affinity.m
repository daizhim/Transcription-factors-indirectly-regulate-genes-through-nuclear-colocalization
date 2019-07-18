clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('TF_binding_ratio.mat');
a=load('test2.txt');
[grow,gcol]=size(a);
[arow,acol]=size(binding2knockout_refined);
k2=0;
k3=0;
ind2=[];
ind3=[];
for i=1:acol
    k=0;
    k1=0;
    ind=[];
    ind1=[];
    [crow,ccol]=size(binding2knockout_refined{i}{1});
    for ii=1:crow
        if (binding2knockout_refined_ind{i}{1}(ii,1)==1 || binding2knockout_refined_ind1{i}{1}(ii,1)==1) && TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1))~=0
            k=k+1;
            ind(k,1)=abs(TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1)));
            k2=k2+1;
            ind2(k2,1)=abs(TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1)));
        elseif  binding2knockout_refined_ind{i}{1}(ii,1)==0 &&  binding2knockout_refined_ind1{i}{1}(ii,1)==0 && TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1))~=0
            k1=k1+1;
            ind1(k1,1)=abs(TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1)));
            k3=k3+1;
            ind3(k3,1)=abs(TF_binding_ratio(binding2knockout_refined{i}{1}(ii,1),a(i,1)));
        end
    end
end
ave1(1,1)=mean(ind2);
bootstat=bootstrp(100,@mean,ind2);
ave1(1,2)=std(bootstat);
ave1(1,3)=mean(ind3);
bootstat=bootstrp(100,@mean,ind3);
ave1(1,4)=std(bootstat);
[p,h]=ranksum(ind2,ind3);
ave1(1,5)=p;
%ave1(1,5) is the resulting P value
