clear;
nobjs = 3;
nweights = 4;
npoints = 10000;
borg_path = 'C:\Users\Barney\Documents\borg\';
ref = 2;
obj_labs = {'J_1','J_2','J_3'};
labs  = {'Random','Borg','Optimal'};

k = nweights - (nobjs - 1);
xopt = rand(nweights,npoints);
xopt(k+1:end,:) = 0.5;

obj_opt = dtlz2(xopt,nobjs)';
hv_opt = Hypervolume_MEX(obj_opt,ones(nobjs,1).*ref);

cols = distinguishable_colors(3);

iter = 20;
obj_rand = cell(iter,1);
obj_borg = cell(iter,1);
hv_rand = zeros(iter,1);
hv_borg = zeros(iter,1);

for i = 1 : iter
    x = rand(nweights,npoints);
    obj = dtlz2(x,nobjs)';
    obj_rand{i} = paretoFront(obj);
    hv_rand(i) = hv_opt - Hypervolume_MEX(obj_rand{i},ones(nobjs,1).*ref);
    obj_borg{i} = dlmread([borg_path num2str(i) '.out']);
    obj_borg{i} = obj_borg{i}(:,nweights+1:end);
    hv_borg(i) = hv_opt - Hypervolume_MEX(obj_borg{i},ones(nobjs,1).*ref);
end
figure;
boxplot([hv_rand,hv_borg,zeros(size(hv_borg))]);
ylim([0,0.06])

figure;
for i = 1 : nobjs
    subplot(1,nobjs,i);
    for j = 1 : iter
        a = sort(obj_borg{j}(:,i));
        plot(a,linspace(0,1,numel(a)),'color',cols(2,:),'linewidth',1); hold on;
        
        a = sort(obj_rand{j}(:,i));
        plot(a,linspace(0,1,numel(a)),'color',cols(1,:),'linewidth',1); hold on;
    end
    
    a = sort(obj_opt(:,i));
    plot(a,linspace(0,1,numel(a)),'color',cols(3,:),'linewidth',2); hold on;
    xlim([0,1.5]);
    axis square
    xlabel(['J' num2str(i)])
    ylabel('CDF');
end
set(gcf,'color','w');

figure;
scatter3(obj_rand{1}(:,1),obj_rand{1}(:,2),obj_rand{1}(:,3),3,cols(1,:),'filled');hold on;
scatter3(obj_borg{1}(:,1),obj_borg{1}(:,2),obj_borg{1}(:,3),3,cols(2,:),'filled');hold on;
scatter3(obj_opt(:,1),obj_opt(:,2),obj_opt(:,3),10,cols(3,:),'filled');hold on;

figure;
combs = [2,1;3,1;3,2];
[b,a] = size(combs);
for i = 1 : b
    subplot(1,b,i);
    for j = 1 : iter
    P = paretoFront(obj_rand{1}(:,combs(i,:)));
    scatter(P(:,1),P(:,2),30,cols(1,:),'filled'); hold on;
    
    P = paretoFront(obj_borg{1}(:,combs(i,:)));
    scatter(P(:,1),P(:,2),30,cols(2,:),'filled'); hold on;
    end
    P = paretoFront(obj_opt(:,combs(i,:)));
    scatter(P(:,1),P(:,2),30,cols(3,:),'filled'); hold on;
end

f = figure;
obj_all = [obj_opt,ones(size(obj_opt,1),1).*3];
obj_all = [obj_all;obj_borg{1},ones(size(obj_borg{1},1),1).*2];
obj_all = [obj_all;obj_rand{1},ones(size(obj_rand{1},1),1).*1];
R = randperm(size(obj_all,1));
for j = 1 : numel(R)
    if obj_all(R(j),nobjs+1) == 3
        lw = 0.05;
    else
        lw = 0.2;
    end
    plot(obj_all(R(j),1:nobjs),'color',cols(obj_all(R(j),nobjs+1),:),'linewidth',lw); hold on;
end
set(f,'visible','off');
set(f, 'Position', [100 100 850 200])
print(f,'myfig','-dpng','-r2000');
