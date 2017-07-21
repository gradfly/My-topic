%∂Øª≠œ‘ æ
figure
filename='wave1D.gif';
h1=plot(x,w_free(1,:),'linewidth',1);
axis([0 L min(min(w_free)) max(max(w_free))])
% axis([0 1 -2 2])
set(gcf,'renderer','zbuffer','position',[423 184 520 300])
set(gca,'xcolor','k','ycolor','k')
set(gcf,'color','w')
grid off
im=frame2im(getframe(gcf));
[A,map]=rgb2ind(im,256);
% pause(0.01)
imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);
for i=1:Ttr
    set(h1,'XData',x,'YData',w_free(i,:));
    title(['t=',num2str(t_tr1(i)),' s'])
    drawnow;pause(tmax/Ttr)%grid off
    im=frame2im(getframe(gcf));
    [A,map]=rgb2ind(im,256);
    imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',tmax/Ttr);
end