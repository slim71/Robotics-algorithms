function animation_sot(hObject, event, robot, pos)
    val = round(get(event.AffectedObject, 'Value'));
    ax = show(robot, pos(:,end,1));
    ax.XLim = [-1,1];
    ax.YLim = [-1,1];
    ax.ZLim = [-1,1];
    disp(['pos= ' num2str(val)])
end