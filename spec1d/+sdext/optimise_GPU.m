function out = optimise_GPU()

g = gpuDevice;
memA = g.AvailableMemory;  % Check available memory
max_length = memA/(8*3);
max_length = max_length/10; % Just to make sure

num_samples = 15;

l = logspace(3,log10(max_length),num_samples);
n = linspace(1E-4,0.25,num_samples);
sdext.setpref('gpuArray',0);

t = zeros(4,num_samples);
% tt = t;

for i = 1:num_samples
    s = spec1d(rand(round(l(i)),1),rand(round(l(i)),1),rand(round(l(i)),1));
    s_g = s.pushToGraphicsCard;
    t(1,i) = gputimeit(@()s_g.combine(n(end),'method','mean'));
    t(3,i) = timeit(@()s.combine(n(end),'method','mean'));
    t(2,i) = gputimeit(@()s_g.rebin(n(end),'method','average'));
    t(4,i) = timeit(@()s.rebin(n(end),'method','average'));
end

% for i = 1:num_samples
%     if i == 1
%         s = spec1d(rand(round(l(end)),1),rand(round(l(end)),1),rand(round(l(end)),1));
%         s_g = s.pushToGraphicsCard;
%     end
%     tt(1,i) = gputimeit(@()s_g.combine(n(i),'method','mean'));
%     tt(3,i) = timeit(@()s.combine(n(i),'method','mean'));
%     tt(2,i) = gputimeit(@()s_g.rebin(n(i),'method','average'));
%     tt(4,i) = timeit(@()s.rebin(n(i),'method','average'));
% end

[~, ind1] = min(abs((t(3,:)./t(1,:) -1)));
[~, ind2] = min(abs((t(4,:)./t(2,:) -1)));

out = round(mean(l([ind1 ind2])));
sdext.setpref('minGPULength',out)

figure;
pl = plot(log10(l),t(3,:)./t(1,:),log10(l),t(4,:)./t(2,:));
legend(pl,{'Combine','Rebin'},'Location','NorthWest')
hold on
plot(log10([out out]),get(gca,'YLim'),'k-')
xlabel('Object size (10^x)')
ylabel('Speedup')

% figure;
% plot(n,tt(3,:)./tt(1,:),n,tt(4,:)./tt(2,:))
% legend({'Combine','Rebin'},'Location','NorthWest')
% xlabel('Bin size')
% ylabel('Speedup')

