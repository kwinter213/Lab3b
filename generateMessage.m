%%Generate Message                      3/12/19
%Kimberly Winter

function ifftMess=generateMessage(preamble, header, message, bufferSize)
    
    totalMessage=[header;message];
    
    ifftMess(bufferSize+1:bufferSize+(64*3))=preamble.';
    for i=1:length(header)/64
        chunk=ifft(totalMessage((i-1)*64+1:i*64));
        ifftMess(bufferSize+(64*3)+(i-1)*80+1:bufferSize+(64*3)+(i-1)*80+16)=chunk(end-15:end);
        ifftMess(bufferSize+(64*3)+(i-1)*80+17:bufferSize+(64*3)+(i*80))=chunk;
    end
    
    for i=length(header)/64+1:length(header)/64+length(message)/64
        chunk=ifft(totalMessage((i-1)*64+1:i*64));
        ifftMess(1000+bufferSize+(64*3)+(i-1)*80+1:1000+bufferSize+(64*3)+(i-1)*80+16)=chunk(end-15:end);
        ifftMess(1000+bufferSize+(64*3)+(i-1)*80+17:1000+bufferSize+(64*3)+(i-1)*80+80)=chunk;
    end
end