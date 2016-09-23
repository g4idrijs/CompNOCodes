function visEchoSubModule(visLocRange, xmt, rcv)

figure
for visLoc = visLocRange'
    % visLoc = position to visualize echo displacement at in meters
    [h,~] = calc_hhp(xmt,rcv, visLoc');   
    
%     if(numel(h) > 100)
%        h = h(1:500); 
%     end
    
    plot(h)
    title(sprintf('Echo at %2.1f mm During First Transmission',visLoc(3)*1000));
    ylabel('Voltage Trace')
    xlabel('Time Index (right is later in time)')
%     ylim([-5.5*1e-23, 5.5*1e-23])
    
    pause(eps)
end
    
end

