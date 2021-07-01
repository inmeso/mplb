
from git.remote import RemoteProgress

class Progress(RemoteProgress):
    def update(self, *args):
        print(self._cur_line)

def UnTar(fileName, dir=""):
    try:
        import tarfile
    except ImportError:
        print("Please install tarfile package first!")
        return None
    print("Extacting the package",fileName)
    if fileName.endswith("tar.gz"):
        tar = tarfile.open(fileName, "r:gz")
    elif fileName.endswith("tar"):
        tar = tarfile.open(fileName, "r:")
    if dir != "":
        tar.extractall(dir)
    else:
        tar.extractall()
    tar.close()

def DownloadUrl(link):
    try:
        import requests
        import urllib
    except ImportError:
        print("Please install requests and urllib package first!")
        return None
    fileName = link.split('/')[-1]
    url = urllib.request.urlopen(link)
    file = open(fileName, 'wb')
    fileSize = int(url.info()["Content-Length"])
    print("Downloading: %s Bytes: %s" % (fileName, fileSize))
    fileSizeDownloaded = 0
    blockSize = 8192
    while True:
        buffer = url.read(blockSize)
        if not buffer:
            break
        fileSizeDownloaded += len(buffer)
        file.write(buffer)
        print('%10d  [%3.2f%%]\r' % (fileSizeDownloaded, fileSizeDownloaded * 100. / fileSize),end="")
    file.close()
    return fileName

def CloneGitRepo(link,dir):
    try:
        from git.repo.base import Repo
    except ImportError:
        print("Please install GitPython package first!")
        return None
    Repo.clone_from(link, dir, progress=Progress())

def RunShellCmd(cmd,silent=False):
    try:
        import subprocess
        import sys
    except ImportError:
        print("Please install subprocess and sys package first!")
        return None
    if isinstance(cmd,str):
        cmd = cmd.split()
    if not silent:
        subprocess.run(cmd, stdout=sys.stdout, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
