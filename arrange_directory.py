#!/usr/bin/su root

import os
import sys
import glob
import pandas as pd
from Bio import Entrez
import subprocess
import mygene
import contextlib2
import numpy as np
import requests, json
from io import StringIO
import datetime
import time


def rename_directory_and_contents(dir_path):
    """
    Rename the specified directory and its contents by removing any part of the name that comes before the first underscore.
    """
    if not os.path.exists(dir_path):
        return "Directory path does not exist."

    # Rename the directory itself
    dir_name = os.path.basename(dir_path)
    underscore_index = dir_name.find('_')
    if underscore_index != -1:
        new_dir_name = dir_name[underscore_index + 1:]
        parent_dir = os.path.dirname(dir_path)
        new_dir_path = os.path.join(parent_dir, new_dir_name)
        os.rename(dir_path, new_dir_path)
        dir_path = new_dir_path  # Update the dir_path to the new path
    else:
        new_dir_name = dir_name

    # Rename files within the directory
    renamed_files = {}
    for file_name in os.listdir(dir_path):
        full_file_path = os.path.join(dir_path, file_name)
        if os.path.isfile(full_file_path):
            underscore_index = file_name.find('_')
            if underscore_index != -1:
                new_file_name = file_name[underscore_index + 1:]
                new_file_path = os.path.join(dir_path, new_file_name)
                os.rename(full_file_path, new_file_path)
                renamed_files[file_name] = new_file_name
            else:
                renamed_files[file_name] = "No change"

    return {"Directory Renamed": new_dir_name, "Files Renamed": renamed_files}


def rename_contents_of_directory(dir_path):
    """
    Rename the contents of the specified directory by removing any part of the name that comes before the first underscore.
    """
    if not os.path.exists(dir_path):
        return "Directory path does not exist."

    renamed_items = {}
    for item_name in os.listdir(dir_path):
        full_item_path = os.path.join(dir_path, item_name)
        underscore_index = item_name.find('_')
        if underscore_index != -1:
            new_item_name = item_name[underscore_index + 1:]
            new_item_path = os.path.join(dir_path, new_item_name)
            os.rename(full_item_path, new_item_path)
            renamed_items[item_name] = new_item_name
        else:
            renamed_items[item_name] = "No change"

    return renamed_items

def main():
    input_dir = "/home/bonus//RNASeq_files/All_OptiDonor_analysis/STAResults/"
    dir_lst = glob.glob(input_dir + "*")

    for dir in dir_lst:
        rename_directory_and_contents(dir)
    dir_lst = glob.glob(input_dir + "*index2")
    dir_lst.remove(input_dir + "BFII112.fastq.gz_index2")
    for dir in dir_lst:
        rename_contents_of_directory(dir)


if __name__ == '__main__':
    main()