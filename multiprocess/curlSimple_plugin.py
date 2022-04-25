#!/usr/bin/env python37

import sys
import argparse
import classad
import io
import json
import os
import pathlib
import posixpath
import pycurl
import shutil
import socket
import time

try:
    from urllib.parse import urlparse # Python 3
except ImportError:
    from urlparse import urlparse # Python 2

DEFAULT_TIMEOUT = 30
CURLSIMPLE_PLUGIN_VERSION = '1.0.0'

def print_help(stream = sys.stderr):
    help_msg = '''Usage: {0} -infile <input-filename> -outfile <output-filename>
       {0} -classad

Options:
  -classad                    Print a ClassAd containing the capablities of this
                              file transfer plugin.
  -infile <input-filename>    Input ClassAd file
  -outfile <output-filename>  Output ClassAd file
  -upload                     Not supported
'''
    stream.write(help_msg.format(sys.argv[0]))

def print_capabilities():
    capabilities = {
         'MultipleFileSupport': False,
         'Version': CURLSIMPLE_PLUGIN_VERSION,
         'PluginType': 'FileTransfer',
         'SupportedMethods': 'curlsimple',
    }
    sys.stdout.write(classad.ClassAd(capabilities).printOld())

def parse_args():
    '''The optparse library can't handle the types of arguments that the file
    transfer plugin sends, the argparse library can't be expected to be
    found on machines running EL 6 (Python 2.6), and a plugin should not
    reach outside the standard library, so the plugin must roll its own argument
    parser. The expected input is very rigid, so this isn't too awful.'''

    # The only argument lists that are acceptable are
    # <this> -classad
    # <this> -infile <input-filename> -outfile <output-filename>
    # <this> -outfile <output-filename> -infile <input-filename>
    if not len(sys.argv) in [2, 5, 6]:
        print_help()
        sys.exit(-1)

    # If -classad, print the capabilities of the plugin and exit early
    if (len(sys.argv) == 2) and (sys.argv[1] == '-classad'):
        print_capabilities()
        sys.exit(0)

    # Upload is not supported
    if '-upload' in sys.argv[1:]:
        print("Upload is not supported. You should not be seeing this.")
        sys.exit(-1)

    # -infile and -outfile must be in the first and third position
    if not (
            ('-infile' in sys.argv[1:]) and
            ('-outfile' in sys.argv[1:]) and
            (sys.argv[1] in ['-infile', '-outfile']) and
            (sys.argv[3] in ['-infile', '-outfile']) and
            (len(sys.argv) == 5)):
        print_help()
        sys.exit(-1)
    infile = None
    outfile = None
    try:
        for i, arg in enumerate(sys.argv):
            if i == 0:
                continue
            elif arg == '-infile':
                infile = sys.argv[i+1]
            elif arg == '-outfile':
                outfile = sys.argv[i+1]
    except IndexError:
        print_help()
        sys.exit(-1)

    return {'infile': infile, 'outfile': outfile}

def format_error(error):
    return '{0}: {1}'.format(type(error).__name__, str(error))

def get_error_dict(error, url = ''):
    error_string = format_error(error)
    error_dict = {
        'TransferSuccess': False,
        'TransferError': error_string,
        'TransferUrl': url,
    }

    return error_dict

def test(debug_type, debug_msg,):
    with open("/home/home4/institut_1b/dmeuser/top_analysis/framework/multiprocess/Output_new.txt", "a") as text_file:
        text_file.write(debug_msg)

class CurlSimplePlugin:

    #def __init__(self):

    def download_data(self, url):

        transfer_stats = []
        
        #  ~curlsimplefind_server = url[(url.find("://") + 3):url.find("?")]
        #  ~curlsimple_args = url[(url.find("?") + 1):].split("&")
        curlsimple_url = url.replace("curlsimple","https")

        # Now transfer each url
        curl = pycurl.Curl()
        start_time = time.time()

        # Setup the transfer
        transfer_success = True
        transfer_error = ""
        file_name = curlsimple_url.split("/")[-1]
        file_size = 0
        file = open(file_name, "wb")
        
        # Use pycurl to actually grab the file
        curl.setopt(curl.URL, curlsimple_url)
        curl.setopt(curl.WRITEDATA, file)
        curl.setopt(curl.CAPATH, os.getenv('X509_CERT_DIR'))
        curl.setopt(curl.SSLCERT, os.getenv('X509_USER_PROXY'))
        curl.setopt(curl.FOLLOWLOCATION, True)
        curl.setopt(curl.VERBOSE, 1)
        #  ~curl.setopt(curl.DEBUGFUNCTION, test)
        try:
            curl.perform()
            file_size = file.tell()
            file.close()
        except Exception as error:
            transfer_success = False
            errno, errstr = error.args
            transfer_error = str(errstr + " (Error " + str(errno) + ")")

        end_time = time.time()

        # Add statistics for this transfer
        this_transfer_stats = {
            'TransferSuccess': transfer_success,
            'TransferProtocol': 'curlsimple',
            'TransferType': 'download',
            'TransferFileName': file_name,
            'TransferFileBytes': file_size,
            'TransferTotalBytes': file_size,
            'TransferStartTime': int(start_time),
            'TransferEndTime': int(end_time),
            'ConnectionTimeSeconds': end_time - start_time,
            'TransferUrl': curlsimple_url,
        }
        if transfer_error:
            this_transfer_stats['TransferError'] = transfer_error

        transfer_stats.append(this_transfer_stats)

        # If this transfer failed, exit immediately
        if not transfer_success:
            curl.close()
            return transfer_stats, False

        # Looks like everything downloaded correctly. Exit success.
        curl.close()
        return transfer_stats, True



if __name__ == '__main__':
    
    debugOut=""
    
    # Start by parsing input arguments
    try:
        args = parse_args()
    except Exception:
        sys.exit(-1)

    curlsimple = CurlSimplePlugin()

    # Parse in the classads stored in the input file. 
    # Each ad represents a single file to be transferred.
    try:
        infile_ads = classad.parseAds(open(args['infile'], 'r'))
    except Exception as err:
        try:
            with open(args['outfile'], 'w') as outfile:
                outfile_dict = get_error_dict(err)
                outfile.write(str(classad.ClassAd(outfile_dict)))
        except Exception:
            pass
        sys.exit(-1)
    
        
    # Now iterate over the list of classads and perform the transfers.
    try:
        with open(args['outfile'], 'w') as outfile:
            for ad in infile_ads:
                try:
                    result_ads, transfer_success = curlsimple.download_data(ad['Url'])
                    # Output the results to file
                    for ad in result_ads:
                        outfile.write(str(classad.ClassAd(ad)))

                    if not transfer_success:
                        sys.exit(-1)

                except Exception as err:
                    print(err)
                    try:
                        outfile_dict = get_error_dict(err, url = ad['Url'])
                        outfile.write(str(classad.ClassAd(outfile_dict)))
                    except Exception:
                        pass
                    sys.exit(-1)

    except Exception:
        sys.exit(-1)
