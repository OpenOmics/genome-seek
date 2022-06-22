# Common helper functions shared across the entire workflow
def provided(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    """

    if not condition:
        # If condition is False, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist


def ignore(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is met (i.e. True) then it will not try to run that rule. This function is the 
    inverse to provided(). 
    """

    if condition:
        # If condition is True, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist


def s3_configured(uri):
    """
    Determines if user can access s3 object using their credentials saved in
    "~/.aws/credentials" or "~/.boto". This handles an edge case where a user
    has aws configure on the target system but the AWS Access Key Id provided in
    their config does not match s3 access policy. This usually occurs when a user
    has setup aws but it is setup for another AWS account or one that does not
    have an IAM role for our S3 bucket.
    :param uri <str>: URI/URL to object in S3 bucket
    :return accessible <boolean>:
        True if user can access S3 object, False if user cannot access the object (403)
    """
    import boto3
    import botocore
    import re

    # Get bucket and key from s3 uri
    parsed = re.match(r's3:\/\/(.+?)\/(.+)', uri)
    bucket, key = parsed.groups()
    accessible = True

    try:
        # Try to access object in s3 bucket
        boto3.resource("s3").Object(bucket, key).load()
    except botocore.exceptions.ClientError as e:
        # User cannot access object in S3 bucket with the credentials
        # stored in "~/.aws/credentials" or "~/.boto".
        if e.response["Error"]["Code"] == "403":
            accessible = False

    return accessible


def abstract_location(file_address, *args, **kwargs):
    """
    Determines if a provided file or list of file(s) resides in a remote location.
    If file(s) are determined to reside in remote store, like a S3 or Google Cloud
    Storage, Snakemake's remote wrapper is used to defined remote files.
    This can be extended further to support more file types listed here.
    https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
    Supported remotes file options include: s3, gs, and sftp
    Input: File path <str> or a list or file paths list[<str>]
    Output: List of files or remote objects
    """

    # Check if user provided any input
    if not file_address or file_address is None:
        raise IOError("Failed to provide any input files! Input(s) are required to resolve required files.".format(file_address))

    # If given file path to one file, convert it a list[<str>]
    file_list = [file_address] if isinstance(file_address, str) else file_address

    # Loop through list of provided files, and if a remote storage option
    # is given, convert its index to a remote file object.
    for i, uri in enumerate(file_list):
        if uri.lower().startswith('s3://'):
            # Remote option for S3 storage
            import snakemake.remote.S3
            import botocore.session

            if botocore.session.get_session().get_credentials() and s3_configured(uri):
                # AWS cli or boto has been configured on target system
                # See ~/.aws/credentials or ~/.boto
                # https://boto.readthedocs.io/en/latest/boto_config_tut.html
                remote_provider = snakemake.remote.S3.RemoteProvider()
            else:
                # If botocore cannot find credentials, try connecting unsigned.
                # This will work for anonymous S3 resources if the resources in the
                # s3 bucket are configured correctly.
                # If a file in provieded as input to a Snakemake rule, only read
                # access is needed to access the remote S3 object.
                remote_provider = snakemake.remote.S3.RemoteProvider(config=botocore.client.Config(signature_version=botocore.UNSIGNED))
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

        elif uri.lower().startswith('gs://'):
            # Remote option for Google Cloud Storage
            import snakemake.remote.GS
            remote_provider = snakemake.remote.GS.RemoteProvider()
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

        elif uri.lower().startswith('sftp://'):
            # Remote option for SFTP transfers
            import snakemake.remote.SFTP
            remote_provider = snakemake.remote.SFTP.RemoteProvider()
            file_list[i] = remote_provider.remote(uri, *args, **kwargs)

    return file_list


def references(config, reflist):
    """
    Checks if a set of required reference files were provided. Some rules depend
    on a set of required reference files that may only exist for specific reference
    genomes. An example of this would be blasklists arriba. The blacklist are manually
    curated and only exist for a few reference genomes (mm10, hg38, hg19).
    If one of the required reference files does not exist, then it will return
    an empty list.
    """

    _all = True
    for ref in reflist:
        try: tmp = config['references'][ref]
        # Check if ref exists in config
        except KeyError:
            _all = False
            break
        # Check if ref is empty key string
        if not tmp: _all = False

    return _all


def allocated(resource, rule, lookup, default="__default__"):
    """Pulls resource information for a given rule. If a rule does not have any information 
    for a given resource type, then it will pull from the default. Information is pulled from
    definitions in the cluster.json (which is used a job submission). This ensures that any 
    resources used at runtime mirror the resources that were allocated.
    :param resource <str>: resource type to look in cluster.json (i.e. threads, mem, time, gres)
    :param rule <str>: rule to lookup its information
    :param lookup <dict>: Lookup containing allocation information (i.e. cluster.json)
    :param default <str>: default information to use if rule information cannot be found
    :return allocation <str>: 
        allocation information for a given resource type for a given rule
    """

    try: 
        # Try to get allocation information
        # for a given rule
        allocation = lookup[rule][resource]
    except KeyError:
        # Use default allocation information
        allocation = lookup[default][resource]
    
    return allocation


def str_bool(s):
    """Converts a string to boolean. It is dangerous to try to
    typecast a string into a boolean value using the built-in 
    `bool()` function. This function avoids any issues that can
    arise when using `bool()`. 
    Example:
      boolean('True') returns True
      boolean('False') returns False
      boolean('asdas') raises TypeError
    """
    val = s.lower()
    if val in ['true', '1', 'y', 'yes']:
        return True
    elif val in ['false', '0', 'n', 'no', '']:
        return False
    else:
        # Provided value could not be
        # type casted into a boolean
        raise TypeError('Fatal: cannot type cast {} into a boolean'.format(val))


def joint_option(prefix, valueslist):
    """Joins a list while adding a common prefix.
    Example:
      joint_option('-i', [1,2,3])
      '-i 1 -i 2 -i 3'
    """
    s = ""
    for v in valueslist:
        s += "{} {} ".format(prefix, v)
    s = s.rstrip()
    return s