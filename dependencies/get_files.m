function new_files = get_files(base_dir, name, excluded_names)
    if nargin < 3; excluded_names = {}; end

    new_files = dir(fullfile(base_dir, name));
    new_files = new_files(~contains({new_files.name}, excluded_names));
    new_files = flip(natsortfiles(new_files));
end
