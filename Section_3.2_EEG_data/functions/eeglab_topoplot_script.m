eeglab;


fileID     = fopen(channels_file,'r');
formatSpec = '%f %f %f';
sizeA      = [3, 129];
A          = fscanf(fileID,formatSpec, sizeA);
fclose(fileID);
A          = A';

fileID = fopen(channels_file_new, 'w');
idx_elec = 1:129;
idx_elec(channels_to_remove) = [];
for i = idx_elec
    fprintf(fileID,'%d\t%f\t%f\t%f\t%s\n', i, -A(i, 2), A(i, 1), A(i, 3), strcat('EEG', num2str(i)));
end
fclose(fileID);

f = figure('visible', 'off');
topoplot(coef_channels, channels_file_new, 'plotrad', plotrad, ...,
    'intrad', intrad, 'headrad', 'rim', 'electrodes', show_electrodes, ...
    'numcontour', numcontour, 'emarker2', ...,
    {cluster_channels,cluster_channels_symbol,cluster_channels_color, ...
    cluster_channels_size});
colormap(color_map);
axis([-0.6 0.6 -0.6 0.6])
colorbar
print(f, '-djpeg', output_channel_figure_file)  
close(f)
