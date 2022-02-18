#!/usr/bin/env python3
# -*- coding: UTF-8 -*-

# Python standard library
from __future__ import print_function
from shutil import copytree
import os, sys, hashlib
import subprocess, json


def md5sum(filename, first_block_only = False, blocksize = 65536):
    """Gets md5checksum of a file in memory-safe manner.
    The file is read in blocks/chunks defined by the blocksize parameter. This is 
    a safer option to reading the entire file into memory if the file is very large.
    @param filename <str>:
        Input file on local filesystem to find md5 checksum
    @param first_block_only <bool>:
        Calculate md5 checksum of the first block/chunk only
    @param blocksize <int>:
        Blocksize of reading N chunks of data to reduce memory profile
    @return hasher.hexdigest() <str>:
        MD5 checksum of the file's contents
    """
    hasher = hashlib.md5()
    with open(filename, 'rb') as fh:
        buf = fh.read(blocksize)
        if first_block_only:
            # Calculate MD5 of first block or chunck of file.
            # This is a useful heuristic for when potentially 
            # calculating an MD5 checksum of thousand or 
            # millions of file.
            hasher.update(buf)
            return hasher.hexdigest()
        while len(buf) > 0:
            # Calculate MD5 checksum of entire file
            hasher.update(buf)
            buf = fh.read(blocksize)

    return hasher.hexdigest()


def permissions(parser, path, *args, **kwargs):
    """Checks permissions using os.access() to see the user is authorized to access
    a file/directory. Checks for existence, readability, writability and executability via:
    os.F_OK (tests existence), os.R_OK (tests read), os.W_OK (tests write), os.X_OK (tests exec).
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        Returns abs path if it exists and permissions are correct
    """
    if not exists(path):
        parser.error("Path '{}' does not exists! Failed to provide vaild input.".format(path))
    if not os.access(path, *args, **kwargs):
        parser.error("Path '{}' exists, but cannot read path due to permissions!".format(path))

    return os.path.abspath(path)


def standard_input(parser, path, *args, **kwargs):
    """Checks for standard input when provided or permissions using permissions().
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param path <str>:
        Name of path to check
    @return path <str>:
        If path exists and user can read from location
    """
    # Checks for standard input
    if not sys.stdin.isatty():
        # Standard input provided, set path as an
        # empty string to prevent searching of '-'
        path = ''
        return path

    # Checks for positional arguments as paths
    path = permissions(parser, path, *args, **kwargs)

    return path


def exists(testpath):
    """Checks if file exists on the local filesystem.
    @param parser <argparse.ArgumentParser() object>:
        argparse parser object
    @param testpath <str>:
        Name of file/directory to check
    @return does_exist <boolean>:
        True when file/directory exists, False when file/directory does not exist
    """
    does_exist = True
    if not os.path.exists(testpath):
        does_exist = False # File or directory does not exist on the filesystem

    return does_exist


def ln(files, outdir):
    """Creates symlinks for files to an output directory.
    @param files list[<str>]:
        List of filenames
    @param outdir <str>:
        Destination or output directory to create symlinks
    """
    # Create symlinks for each file in the output directory
    for file in files:
        ln = os.path.join(outdir, os.path.basename(file))
        if not exists(ln):
                os.symlink(os.path.abspath(os.path.realpath(file)), ln)


def which(cmd, path=None):
    """Checks if an executable is in $PATH
    @param cmd <str>:
        Name of executable to check
    @param path <list>:
        Optional list of PATHs to check [default: $PATH]
    @return <boolean>:
        True if exe in PATH, False if not in PATH
    """
    if path is None:
        path = os.environ["PATH"].split(os.pathsep)

    for prefix in path:
        filename = os.path.join(prefix, cmd)
        executable = os.access(filename, os.X_OK)
        is_not_directory = os.path.isfile(filename)
        if executable and is_not_directory:
            return True
    return False


def err(*message, **kwargs):
    """Prints any provided args to standard error.
    kwargs can be provided to modify print functions 
    behavior.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    print(*message, file=sys.stderr, **kwargs)



def fatal(*message, **kwargs):
    """Prints any provided args to standard error
    and exits with an exit code of 1.
    @param message <any>:
        Values printed to standard error
    @params kwargs <print()>
        Key words to modify print function behavior
    """
    err(*message, **kwargs)
    sys.exit(1)


def require(cmds, suggestions, path=None):
    """Enforces an executable is in $PATH
    @param cmds list[<str>]:
        List of executable names to check
    @param suggestions list[<str>]:
        Name of module to suggest loading for a given index
        in param cmd.
    @param path list[<str>]]:
        Optional list of PATHs to check [default: $PATH]
    """
    error = False
    for i in range(len(cmds)):
        available = which(cmds[i])
        if not available:
            error = True
            err("""\x1b[6;37;41m\n\tFatal: {} is not in $PATH and is required during runtime!
            └── Solution: please 'module load {}' and run again!\x1b[0m""".format(cmds[i], suggestions[i])
            )

    if error: fatal()

    return 


def safe_copy(source, target, resources = []):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            copytree(os.path.join(source, resource), destination)


def git_commit_hash(repo_path):
    """Gets the git commit hash of the RNA-seek repo.
    @param repo_path <str>:
        Path to RNA-seek git repo
    @return githash <str>:
        Latest git commit hash
    """
    githash = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd = repo_path).strip()
    # Typecast to fix python3 TypeError (Object of type bytes is not JSON serializable)
    # subprocess.check_output() returns a byte string
    githash = str(githash)

    return githash


def join_jsons(templates):
    """Joins multiple JSON files to into one data structure
    Used to join multiple template JSON files to create a global config dictionary.
    @params templates <list[str]>:
        List of template JSON files to join together
    @return aggregated <dict>:
        Dictionary containing the contents of all the input JSON files
    """
    # Get absolute PATH to templates in rna-seek git repo
    repo_path = os.path.dirname(os.path.abspath(__file__))
    aggregated = {}

    for file in templates:
        with open(os.path.join(repo_path, file), 'r') as fh:
            aggregated.update(json.load(fh))

    return aggregated


def check_cache(parser, cache, *args, **kwargs):
    """Check if provided SINGULARITY_CACHE is valid. Singularity caches cannot be
    shared across users (and must be owned by the user). Singularity strictly enforces
    0700 user permission on on the cache directory and will return a non-zero exitcode.
    @param parser <argparse.ArgumentParser() object>:
        Argparse parser object
    @param cache <str>:
        Singularity cache directory
    @return cache <str>:
        If singularity cache dir is valid
    """
    if not exists(cache):
        # Cache directory does not exist on filesystem
        os.makedirs(cache)
    elif os.path.isfile(cache):
        # Cache directory exists as file, raise error
        parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
        The provided --singularity-cache already exists on the filesystem as a file.
        Please run {} again with a different --singularity-cache location.
        """.format(sys.argv[0]))
    elif os.path.isdir(cache):
        # Provide cache exists as directory
        # Check that the user owns the child cache directory
        # May revert to os.getuid() if user id is not sufficent
        if exists(os.path.join(cache, 'cache')) and os.stat(os.path.join(cache, 'cache')).st_uid != os.getuid():
                # User does NOT own the cache directory, raise error
                parser.error("""\n\t\x1b[6;37;41mFatal: Failed to provided a valid singularity cache!\x1b[0m
                The provided --singularity-cache already exists on the filesystem with a different owner.
                Singularity strictly enforces that the cache directory is not shared across users.
                Please run {} again with a different --singularity-cache location.
                """.format(sys.argv[0]))

    return cache


def unpacked(nested_dict):
    """Generator to recursively retrieves all values in a nested dictionary.
    @param nested_dict dict[<any>]:
        Nested dictionary to unpack
    @yields value in dictionary 
    """
    # Iterate over all values of given dictionary
    for value in nested_dict.values():
        # Check if value is of dict type
        if isinstance(value, dict):
            # If value is dict then iterate over
            # all its values recursively
            for v in unpacked(value):
                yield v
        else:
            # If value is not dict type then
            # yield the value
            yield value


if __name__ == '__main__':
    # Calculate MD5 checksum of entire file 
    print('{}  {}'.format(md5sum(sys.argv[0]), sys.argv[0]))
    # Calcualte MD5 cehcksum of 512 byte chunck of file,
    # which is similar to following unix command: 
    # dd if=utils.py bs=512 count=1 2>/dev/null | md5sum 
    print('{}  {}'.format(md5sum(sys.argv[0], first_block_only = True, blocksize = 512), sys.argv[0]))
