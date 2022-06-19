function [A] = StateEq(tk,tl)

dt = tk-tl;
A = [1,dt,dt^2/2,0;0,1,dt,0;0,0,1,0;0,0,0,1];

end