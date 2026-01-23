"""Save files from urls.

Using the url to download/save files in .csv, .txt, .jpg, .png, .h5, .mat formats.

Example using requests:

>>> import requests
>>> print('Beginning file download with requests')
>>> url = 'http://i3.ytimg.com/vi/J---aiyznGQ/mqdefault.jpg'
>>> r = requests.get(url)
>>> with open('cat3.jpg', 'wb') as f:
>>>     f.write(r.content)
>>> # Retrieve HTTP meta-data
>>> print(r.status_code)
>>> print(r.headers['content-type'])
>>> print(r.encoding)

Example using urllib.request:

>>> import urllib.request
>>> urllib.request.urlretrieve("https://docs.scipy.org/doc/scipy/reference/_images/interpolate-3.png", "filename.png")

Example 3:

>>> import os
>>> import requests
>>> url = 'https://apod.nasa.gov/apod/image/1701/potw1636aN159_HST_2048.jpg'
>>> page = requests.get(url)
>>> f_ext = os.path.splitext(url)[-1]
>>> f_name = 'img{}'.format(f_ext)
>>> with open(f_name, 'wb') as f:
>>>     f.write(page.content)

"""

import os
import urllib.request


def save_file(url, filename=None):
    """save file from url

    Download/save file from web to local.

    Args:
        url (str): the url for the file you want to save.
        filename (str): optional, use 'filename' if not specified.

    Returns:
        file saved in specified directory.

    """
    if filename:
        urllib.request.urlretrieve(
            url, os.path.join(os.path.dirname(__file__), filename + url.split(".")[-1])
        )
    else:
        urllib.request.urlretrieve(
            url, os.path.join(os.path.dirname(__file__), url.split("/")[-1])
        )


# example:
url = "https://docs.scipy.org/doc/scipy/reference/_images/interpolate-3.png"

filename = url.split("/")[-1]
# print(filename)

datadir = os.path.join(os.path.dirname(__file__), "data", filename)

urllib.request.urlretrieve(url, datadir)

#######################################################################
url1 = "https://upload.wikimedia.org/wikipedia/commons/2/20/Google-Logo.svg"
resource = urllib.request.urlopen(url1)
output = open(
    os.path.join(os.path.dirname(__file__), "data", url1.split("/")[-1]), "wb"
)
output.write(resource.read())
output.close()

url2 = "https://www.fu-berlin.de/sites/corporate-design/grundlagen/_medien/fu_logo.png"
resource = urllib.request.urlopen(url2)
output = open(
    os.path.join(os.path.dirname(__file__), "data", url2.split("/")[-1]), "wb"
)
output.write(resource.read())
output.close()
#######################################################################


"""
def form_filename(i):
    length_i = len(i)
    if length_i == 1:
        return '000' + str(i)
    if length_i == 2:
        return '00' + str(i)
    if length_i == 3:
        return '0' + str(i)

M = 1
N = 10
for i in range(M, N):
    filename = form_filename(str(i)) + ".jpg"
    # filename = '5544.jpg'
    url = "https://i.guim.co.uk/img/media/7a633730f5f90db3c12f6efc954a2d5b475c3d4a/0_138_5544_3327/master/" + filename
    urllib.request.urlretrieve(url, filename)

"""
