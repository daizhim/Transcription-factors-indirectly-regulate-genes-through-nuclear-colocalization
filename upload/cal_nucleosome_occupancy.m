clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('transcript2nucleosome_lee.mat');
load('TF_noboundgene_MacIsaac.mat');
load('TSS.mat');
[arow,acol]=size(binding2knockout_refined);
[drow,dcol]=size(transcript2nucleosome);
val=zeros(6767,2);
val(:,1)=1:6767;
val(:,2)=-10000;
for i=1:dcol
    if nnz(transcript2nucleosome{i}.utr5(1:500,1))~=0
        val(TSS(i,5),2)=sum(transcript2nucleosome{i}.utr5(1:500,1))/nnz(transcript2nucleosome{i}.utr5(1:500,1));
    end
end
k=0;
k1=0;
ind=[];
ind1=[];
for i=1:acol
    [crow,ccol]=size(binding2knockout_refined{i}{2});
    for ii=1:crow
        if binding2knockout_refined_ind{i}{2}(ii,1)==1 || binding2knockout_refined_ind1{i}{2}(ii,1)==1
            k=k+1;
            ind(k,1)=binding2knockout_refined{i}{2}(ii,1);
        end
    end
    [crow,ccol]=size(binding2knockout_refined{i}{1});
    ind1(k1+1:k1+crow,1)=binding2knockout_refined{i}{1};
    k1=k1+crow;
end
ind=sortrows(ind,1);
new_ind=[];
new_ind(1,1)=ind(1,1);
[drow,dcol]=size(ind);
k=1;
for i=2:drow
    if ind(i,1)~=new_ind(k,1)
        k=k+1;
        new_ind(k,1)=ind(i,1);
    end
end
ind=new_ind;

ind1=sortrows(ind1,1);
new_ind=[];
new_ind(1,1)=ind1(1,1);
[drow,dcol]=size(ind1);
k=1;
for i=2:drow
    if ind1(i,1)~=new_ind(k,1)
        k=k+1;
        new_ind(k,1)=ind1(i,1);
    end
end
ind1=new_ind;

[drow,dcol]=size(ind);
[erow,ecol]=size(ind1);
k=0;
ind2=[];
for i=1:drow
    temp_ind=0;
    for j=1:erow
        if ind(i,1)==ind1(j,1)
            temp_ind=1;
            break;
        end
    end
    if temp_ind==0
        k=k+1;
        ind2(k,1)=ind(i,1);
    end
end

[drow,dcol]=size(ind);
[erow,ecol]=size(ind1);
k=0;
ind3=[];
for i=1:erow
    temp_ind=0;
    for j=1:drow
        if ind1(i,1)==ind(j,1)
            temp_ind=1;
            break;
        end
    end
    if temp_ind==0
        k=k+1;
        ind3(k,1)=ind1(i,1);
    end
end

[drow,dcol]=size(ind2);
k=0;
ind4=[];
for i=1:6767
    temp_ind=0;
    for j=1:drow
        if i==ind2(j,1)
            temp_ind=1;
            break;
        end
    end
    if temp_ind==0
        k=k+1;
        ind4(k,1)=i;
    end
end
k=0;
value=[];
[drow,dcol]=size(ind2);
for j=1:drow
    if val(ind2(j,1),2)~=-10000
        k=k+1;
        value(k,1)=val(ind2(j,1),2);
    end
end
k=0;
value1=[];
[drow,dcol]=size(ind3);
for j=1:drow
    if val(ind3(j,1),2)~=-10000
        k=k+1;
        value1(k,1)=val(ind3(j,1),2);
    end
end
k=0;
value2=[];
[drow,dcol]=size(ind4);
for j=1:drow
    if val(ind4(j,1),2)~=-10000
        k=k+1;
        value2(k,1)=val(ind4(j,1),2);
    end
end
ave(1,1)=mean(value);
bootstat=bootstrp(100,@mean,value);
ave(1,2)=std(bootstat);
ave(1,3)=mean(value1);
bootstat=bootstrp(100,@mean,value1);
ave(1,4)=std(bootstat);
ave(1,5)=mean(value2);
bootstat=bootstrp(100,@mean,value2);
ave(1,6)=std(bootstat);
[p,h]=ranksum(value,value1);
pvalue(1,1)=p;
[p,h]=ranksum(value,value2);
pvalue(1,2)=p;
%pvalue is the resulting P value
