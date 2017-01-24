function save_figure_to_file(figure_handle, figure_save_filepath, figure_name)

% This function saves a figure to file in .fig, .eps and .png formats.

% INPUTS:
% figure_handle: Object storing the figure properties. Obtained by calling gcf.
% figure_save_filepath: Path to the folder where the figure will be saved
% figure_name: The name under which the figure will be saved

% OUTPUTS:
% This function does not return any output

% AUTHOR:
% Lalit Rajendran (lrajendr@purdue.edu)

savefig(figure_handle, [figure_save_filepath figure_name '.fig']);
print(figure_handle, [figure_save_filepath figure_name '.eps'], '-depsc');
print(figure_handle, [figure_save_filepath figure_name '.png'], '-dpng');

end