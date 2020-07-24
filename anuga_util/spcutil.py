"""
Create output directories if not exists
"""
from os import path, mkdir, sep


def dirStruct(prj_name, main_output_dir='.' + sep + 'output'):
    if not path.exists(main_output_dir):
        mkdir(main_output_dir)
        print('Main output directory: ', main_output_dir, " created under project directory")
    else:
        print('Main output directory: ', main_output_dir, " exists")

    if not path.exists(main_output_dir + sep + prj_name):
        mkdir(main_output_dir + sep + prj_name)
        print(prj_name, " created under ", main_output_dir, " directory")
    else:
        raise Exception('Name conflict. Change suffix-prefix or delete/move ' + prj_name)

    if not path.exists(main_output_dir + sep + prj_name + sep + 'gauges'):
        mkdir(main_output_dir + sep + prj_name + sep + 'gauges')
        print(prj_name + sep + 'gauges', " created under ", main_output_dir, " directory")
