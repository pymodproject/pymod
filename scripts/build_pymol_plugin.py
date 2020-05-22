import os
import argparse
import datetime
import zipfile


# Parses the commandline.
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--revision", help="Number of the plugin revision for this day.", type=int, default=None)
parser.add_argument("-p", "--pymod_plugin_dirpath", help="Path of th 'pymod3'.", type=str, default="pymod3")
cmd = parser.parse_args()


# Name of the PyMOL plugin file .zip file.
if cmd.revision is None:
    plugin_archive_name = "pymod3.zip"
else:
    current_date = datetime.datetime.now()
    plugin_archive_name = "pymod3_%s_%s_%s_rev_%s.zip" % (current_date.day,
                                                          current_date.month,
                                                          current_date.year,
                                                          cmd.revision)
target_filepath = plugin_archive_name


print("- Building a PyMOL plugin file from directory: %s" % cmd.pymod_plugin_dirpath)
walk_data = os.walk(cmd.pymod_plugin_dirpath)

with zipfile.ZipFile(target_filepath, 'w', zipfile.ZIP_DEFLATED) as z_fh:
    for root, dirs, files in walk_data:
        for filename in files:
            if filename.endswith(".pyc"):
                continue
            z_fh.write(os.path.join(root, filename))

print("- Completed. PyMOL plugin file written: %s" % target_filepath)
