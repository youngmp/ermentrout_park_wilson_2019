function [value,isterminal,direction] = SpikeEvent(t,y)

value = y(1)-30;
isterminal=1;
direction=1;

end