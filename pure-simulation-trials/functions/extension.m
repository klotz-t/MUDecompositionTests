function eY = extension(Y,extfact)
% Extend a multi-channel signal 'Y' by an extansiofactor 'extfact'.
eY = zeros(size(Y,1)*extfact,size(Y,2)+extfact);
for index = 1:extfact
    eY((index-1)*size(Y,1)+1:index*size(Y,1),[1:size(Y,2)]+(index-1)) = Y;
end
end