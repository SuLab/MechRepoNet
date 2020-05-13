import os
import sys
import ftplib
import urllib
import argparse
import requests
import pandas as pd

def get_content_type(url):
    """
    Determines the content type of a URL
    """
    h = requests.head(url, allow_redirects=True)
    header = h.headers
    content_type = header.get('content-type')
    return content_type


def is_downloadable(url):
    """
    Determines if url is a downloadable resource
    """
    content_type = get_content_type(url)
    if is_text(url, content_type):
        return False
    if 'html' in content_type.lower():
        return False
    return True


def is_text(url, content_type=None):
    """
    Determines if a URL is text
    """
    # Allow user to skip step if content type already retrevied
    if content_type is None:
        content_type = get_content_type(url)
    # some text files have no content type header
    if content_type is None or "text" in content_type:
        return True
    return False


def download_file(url, out_name):
    """
    Downloads a file ot the given directory
    """
    r = requests.get(url, allow_redirects=True)
    open(out_name, 'wb').write(r.content)


def save_text(url, out_name):
    """Save text to given filename"""
    r = requests.get(url)
    open(out_name, 'w').write(r.text)


def parse_path(path):
    """Parse a unix path with relative and/or home directory markers"""
    out_path = path
    if '~' in out_path:
        out_path = os.path.expanduser(out_path)
    if '.' in out_path:
        out_path = os.path.realpath(os.path.abspath(out_path))
    else:
        out_path = os.path.abspath(out_path)
    return out_path


def is_ftp(url):
    """Determine if a url is over ftp protocol"""
    return urllib.parse.urlparse(url).scheme == 'ftp'


def download_ftp(url, out_name):
    """Download a file from an ftp server"""
    # Parse the FTP url
    parsed = urllib.parse.urlparse(url)
    server = parsed.netloc
    dl_path = parsed.path

    # Download the file
    ftp = ftplib.FTP(server)
    ftp.login()
    ftp.retrbinary("RETR {}".format(dl_path), open(out_name, 'wb').write)
    ftp.quit()


if __name__ == "__main__":
    # Get command-line arguments
    parser = argparse.ArgumentParser(description='Download Files Needed for HetNet Building')
    parser.add_argument('-r', '--re_download', help='Re-Download files that already exist in the out locaiton',
                        action='store_true')
    parser.add_argument('-o', '--out_dir', help='Specifiy a non-default directory to save files to', type=str,
                        default='')

    # Parse the agruments
    args = parser.parse_args()
    re_download = args.re_download
    out_dir = args.out_dir

    # Read file directory
    read_dir = os.path.abspath('../0_data/manual')

    # Get default output directory if none provided
    if not out_dir:
        f_name = sys.argv[0]
        base_name = os.path.splitext(f_name)[0]
        out_dir = os.path.join(os.path.abspath('../2_pipeline'), base_name, 'out')

    # Make the write directory if it does not already exist
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)

    # Read in info on files to download
    file_data = pd.read_csv(os.path.join(read_dir, 'edge_file_info.csv'))

    # Download the files
    for row in file_data.itertuples():
        out_name = os.path.join(out_dir, row.file_name)

        if os.path.exists(out_name) and not re_download:
            print('File {} exits. Skipping...'.format(row.file_name))
            continue

        if is_ftp(row.url):
            print('Getting {} from ftp server'.format(row.file_name))
            download_ftp(row.url, out_name)
            print('Done')
        elif is_downloadable(row.url):
            print('Downloading {}'.format(row.file_name))
            download_file(row.url, out_name)
            print('Done')
        elif is_text(row.url):
            print('Saving {}'.format(row.file_name))
            save_text(row.url, out_name)
            print('Done')
        else:
            print(row.file_name, ": Not a downloadable file or text... ")
            print('Skipping....')

