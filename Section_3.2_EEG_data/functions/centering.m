function [Xc] = centering(X,mode)
%input: tensor X of size ijk, and mode [a,b or c] default c (ixj)
%output: centered tensor X

[i,j,k] = size(X);

if nargin<2,
   mode='c';
end
    
switch mode % mode is a string
    
    case{'a'} % through mode a (jxk) :
        
        for a=1:i
            x=X(a,:,:);
            x=reshape(x,j,k);
            [m,n]=size(x);
            mx=mean(x);
            xc=(x-mx(ones(m,1),:));
            Xc(i,:,:)=    xc(:,:);
        end
        
    case{'b'}  % through mode a (ixk) :
        
        for b=1:j
            x=X(:,b,:);
            x=reshape(x,i,k);
            [m,n]=size(x);
            mx=mean(x);
            xc=(x-mx(ones(m,1),:));
            Xc(:,b,:)=    xc(:,:);
        end
        
    case{'c'} % through mode a (ixj):
        
        for c=1:k
            x=X(:,:,c);
            [m,n]=size(x);
            mx=mean(x);
            xc=(x-mx(ones(m,1),:));
            Xc(:,:,c)=    xc(:,:);
        end
        
    otherwise
        error('please choose between mode a,b or c')
end
