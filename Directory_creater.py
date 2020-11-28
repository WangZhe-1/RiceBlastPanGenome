from pathlib import Path
def directory_creater(directory_name):
    directory_path=Path(directory_name)
    if directory_path.is_dir():
        return((directory_path))
    else:
        directory_path.mkdir()
        return(directory_path)