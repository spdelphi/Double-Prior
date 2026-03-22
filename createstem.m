function figure1 = createstem(x,abundance,tle)
%CREATEFIGURE(ymatrix1)
%  YMATRIX1:  stem 矩阵数据

%  由 MATLAB 于 07-Dec-2024 18:21:42 自动生成
ymatrix1 = [x;abundance]';

% 创建 figure
figure1 = figure;

% 创建 axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% 使用 stem 的矩阵输入创建多行
stem1 = stem(ymatrix1);

% h_legend = legend('');
% set(h_legend, 'String', {'sum(|\hat X|,2)','sum(|X_E|,2)'})

set(stem1(1),'DisplayName','sum(|X_T|,2)');
set(stem1(2),'DisplayName','sum(|X_E|,2)');
ylim([0,500]);
% 创建 title
title({tle},'FontWeight','bold','FontSize',16,'FontName','Times New Roman');
xlabel('Endmembers','FontSize',16,'FontName','Times New Roman');
box(axes1,'on');
hold(axes1,'off');
% 创建 legend
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.611904761904763 0.798809523809523 0.201785714285714 0.112301587301587]);

