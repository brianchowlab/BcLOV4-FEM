load('vcell_TIRF_result.mat')
t = 0:0.5:100;
plot(t,max_height_vcell,'k','LineWidth',2)
f.Position = [1000 500 700 600];
ax = gca;
ax.TickLength = [0.025, 0.025];
xlabel('Time (s)')
ylabel('Amplitude (molecules/um^2)')
set(gca,'FontSize',30)
box on
set(gca,'linew',2)
set(gca,'tickdir','out')
box off
xlim([0,60])
