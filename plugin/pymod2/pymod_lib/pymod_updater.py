import os
import sys
import shutil
import zipfile
import urllib
import pymod_os_specific as pmos


url_base = "http://schubert.bio.uniroma1.it/downloads/"
plugin_zipfile_name = "pymod_stable.zip"
plugin_releasefile_name = "pymod_stable.txt"
print_info = False
# urllib.urlretrieve(url_base+"/"+code+".pdb", code+".pdb")


def check_for_updates(pymod_version, pymod_release):

    if print_info:
        print "- Checking for PyMod updates."

    # Gets the release of the PyMod version currently in use.
    local_pymod_release = pymod_version+"."+pymod_release
    local_pymod_release_float = local_pymod_release.split(".")[0]+"."+"".join(local_pymod_release.split(".")[1:])

    # Gets the 'pymod_latest.txt' file which contains information about the latest PyMod release
    # number.
    plugin_releasefile_temp_name = urllib.urlretrieve(url_base+plugin_releasefile_name)[0]
    vfh = open(plugin_releasefile_temp_name,"r")
    network_pymod_release = vfh.read().rstrip()
    network_pymod_release_float = network_pymod_release.split(".")[0]+"."+"".join(network_pymod_release.split(".")[1:])
    vfh.close()

    if print_info:
        print "- Latest PyMod release from network:", network_pymod_release_float
        print "- Currently installed PyMod release:", local_pymod_release_float

    if network_pymod_release_float > local_pymod_release_float:
        return network_pymod_release
    else:
        return False


def fetch_plugin_zipfile():
    if print_info:
        print "- Fetching"
    # Gets the PyMod plugin zip file latest stable release.
    plugin_zipfile_temp_name = urllib.urlretrieve(url_base+plugin_zipfile_name, "temp_plugin_zipfile.zip")[0]
    return plugin_zipfile_temp_name


def update_pymod(plugin_zipfile_temp_name, pymod_plugin_dir):
    """
    Called in 'pymod_main.py' in order to update the plugin files to the latest stable version.
    """
    try:
        # Unpacks the zipfile in a temp directory.
        stable_pymod_dir_temp_name = "stable_pymod_dir_temp"
        zfh = open(plugin_zipfile_temp_name, 'rb')
        zipfile_obj = zipfile.ZipFile(zfh)
        pmos.zipfile_extract_all(zipfile_obj, stable_pymod_dir_temp_name)
        zfh.close()

        # Copies the new files to the PyMod plugin directory by overwriting the old ones.
        target_dir = pymod_plugin_dir
        for path, dirs, files in os.walk(stable_pymod_dir_temp_name):
            rpath = os.path.sep.join(path.split(os.path.sep)[1:])
            for d in dirs:
                if not os.path.isdir(os.path.join(target_dir,rpath, d)):
                    os.mkdir(os.path.join(target_dir,rpath, d))
            for f in files:
                if os.path.isfile(os.path.join(target_dir,rpath, f)):
                    os.remove(os.path.join(target_dir,rpath, f))
                shutil.move(os.path.join(path,f), os.path.join(target_dir,rpath, f))
        return (True, "Update Successful")

    except Exception, e:
        # Remove temp files.
        if os.path.isfile(plugin_zipfile_temp_name):
            os.remove(plugin_zipfile_temp_name)
        if os.path.isdir(stable_pymod_dir_temp_name):
            shutil.rmtree(stable_pymod_dir_temp_name)
        return (False, "PyMod update failed because of the following error: %s" % e)
