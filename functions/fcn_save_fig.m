function full_filename_with_path=fcn_save_fig(file_name_prefix,save_folder,fig_file_extension,overwrite_flag)

full_filename_with_path=strcat(save_folder,file_name_prefix,fig_file_extension);

if isempty(overwrite_flag)
    if exist(full_filename_with_path,'file')==0
        export_fig(full_filename_with_path,'-transparent','-nocrop')
    end

else
        export_fig(fig_name,'-transparent','-nocrop')
end
