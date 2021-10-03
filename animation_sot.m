function animation_sot(hObject, event, robot, pos)
    val = round(get(event.AffectedObject, 'Value'));
    ax = show(robot, pos(:,end,1));
    ax.XLim = [0, 0.8];%[-1,1];
    ax.YLim = [0, 0.8];%[-1,1];
    ax.ZLim = [0, 1]; %[-1,1];
    disp(['pos= ' num2str(val)])
end