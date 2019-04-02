%%Kimberly Winter                       4/2/19
%Main function OFDM with Schmidl-cox

%Generate message to send through channel using generateMessage
bufferSize=50;

header=generateRand(64*100);
message=generateRand(19200);
preamble=generateRand(64);
preamble=[preamble; preamble; preamble];
mess2send=generateMessage(preamble, header,message,bufferSize);

%Send message through channel

%plot(real(mess2send));
receivedMess=nonflat_channel_timing_error(mess2send.');

%Estimate H
%Trim header and divide by known header
[correlation,lag] = xcorr(receivedMess,preamble.');
[M,I]=max(abs(correlation));
lagDiff=lag(I);
receivedMess=receivedMess.';

%trim preamble
trimmedPreamble=receivedMess(lagDiff+1:lagDiff+(64*3));

%estimate fdelta for each symbol in preamble
est=zeros(64,1);
for i=64:127
    est(i-63)=angle(trimmedPreamble(i+64)/trimmedPreamble(i));
end

%take average of all fdeltas
fdelta_est=mean(est)/64

%Adjust for fdelta
for i=1:length(receivedMess)
   adjustedMess(i)=receivedMess(i)/exp(j*i*fdelta_est); 
end

%Change indices here
trimmedHead=adjustedMess(lagDiff+(64*3)+1:lagDiff+(64*3)+length(header)*80/64);
trimmedMessage=adjustedMess(lagDiff+(64*3)+length(header)*80/64+1000+1:lagDiff+(64*3)+length(header)*80/64+1000+length(message)*80/64);

%take FFT of trimmed signal for header
for i=1:(length(trimmedHead))/80
    final(64*(i-1)+1:64*(i-1)+64)= fft(trimmedHead((i-1)*80+17:(i-1)*80+80));
end

for i=1:(length(trimmedMessage))/80
    final(6400+1+64*(i-1):6400+64*(i-1)+64)= fft(trimmedMessage((i-1)*80+17:(i-1)*80+80));
end

%Divide Yk/Hk
HEstimate=mean(final(1:6400)./header.');

MessageEstimate=final(6401:6400+length(message))/HEstimate;
normalEst=normalize(MessageEstimate);

error=0;
for i=1:length(normalEst)
    if(normalEst(i)~=message(i))
        error=error+1;
        i
    end
end