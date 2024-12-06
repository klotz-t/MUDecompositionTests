function eY = extension(Y,R)
% Extend a multi-channel signal 'Y' by an extansiofactor 'R'.
% eY = zeros(size(Y,1)*R,size(Y,2)+R);
% for index = 1:R
%     eY((index-1)*size(Y,1)+1:index*size(Y,1),[1:size(Y,2)]+(index-1)) = Y;
% end

[n,N] = size(Y);
eY = zeros(n*R,N);
for i=1:n
    eY((i-1)*R+1:i*R,:)=toeplitz([Y(i,1);zeros(R-1,1)],Y(i,:));
end

end