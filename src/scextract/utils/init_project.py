import os

def init_project(config_path='config.ini',
                 overwrite=False):
                 
    # Check if config file exists
    if os.path.exists(config_path) and not overwrite:
        raise ValueError('Config file already exists. Set overwrite=True to overwrite.')

    # Create config file from template
    os.system(f'cp {os.path.dirname(__file__)}/config_template.ini {config_path}')
    print(f'Config file created at {config_path}')