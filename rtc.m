function [ w2 ] = rtc( w,x,y )
[m,n]=size(w);
w2=zeros(m,n);
for i=1:m
    w_temp=reshape(w(i,:),x,y);
    w2(i,:)=reshape(w_temp',1,x*y);
end

