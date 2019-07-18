clear
load('binding2knockout_Reimand3.mat');
binding2knockout_refined_ind1=binding2knockout_refined_ind;
load('binding2knockout_Reimand3_inter.mat');
load('YPD_correlation1.mat');
[arow,acol]=size(result);
[brow,bcol]=size(expression_cor);
temp_expression_cor1=zeros(brow*brow,1);
temp_expression_cor2=temp_expression_cor1;
temp_expression_cor3=temp_expression_cor1;
temp_expression_cor4=temp_expression_cor1;
k=0;
for i=1:brow
    for j=i+1:bcol
        if expression_cor(i,j)~=0
            k=k+1;
            temp_expression_cor1(k,1)=expression_cor(i,j);
        end
    end
end
expression_cor1=temp_expression_cor1(1:k,1);
ind1=0;
ind2=0;
ind3=0;
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
        else
            k2=k2+1;
            temp_ind2(k2,1)=binding2knockout_refined{i}{2}(j,1);
        end
    end
    [crow,ccol]=size(binding2knockout{i}{1});
    k4=0;
    temp_cor3=[];
    for ii=1:k1
        for jj=1:crow
            if expression_cor(temp_ind1(ii,1),binding2knockout{i}{1}(jj,1))~=0
                k4=k4+1;
                temp_cor3(k4,1)=expression_cor(temp_ind1(ii,1),binding2knockout{i}{1}(jj,1));
                ind2=ind2+1;
                temp_expression_cor3(ind2,1)=expression_cor(temp_ind1(ii,1),binding2knockout{i}{1}(jj,1));
            end
        end
    end
    k5=0;
    temp_cor4=[];
    for ii=1:k2
        for jj=1:crow
            if expression_cor(temp_ind2(ii,1),binding2knockout{i}{1}(jj,1))~=0
                k5=k5+1;
                temp_cor4(k5,1)=expression_cor(temp_ind2(ii,1),binding2knockout{i}{1}(jj,1));
                ind3=ind3+1;
                temp_expression_cor4(ind3,1)=expression_cor(temp_ind2(ii,1),binding2knockout{i}{1}(jj,1));
            end
        end
    end
end
expression_cor3=temp_expression_cor3(1:ind2,1);                
expression_cor4=temp_expression_cor4(1:ind3,1);   
ave1(1,1)=mean(expression_cor3);
bootstat=bootstrp(100,@mean,expression_cor3);
ave1(1,2)=std(bootstat);
ave1(1,3)=mean(expression_cor4);
bootstat=bootstrp(100,@mean,expression_cor4);
ave1(1,4)=std(bootstat);
[pvalue,h]=ranksum(expression_cor4,expression_cor3);
%pvalue is the resulting P value