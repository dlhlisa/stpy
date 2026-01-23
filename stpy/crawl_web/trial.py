import os
import urllib.request

# example:
# figure out about the path for file saving
print("basename: ", os.path.basename(__file__))
print("dirname: ", os.path.dirname(__file__))
print("abspath: ", os.path.abspath(__file__))
print("abs dirname: ", os.path.dirname(os.path.abspath(__file__)))
print("getcwd: ", os.getcwd())

"""
>>> print('[change directory]')
>>> os.chdir(os.path.dirname(os.path.abspath(__file__)))
>>> print('getcwd: ', os.getcwd())
"""

# example:
url = "https://docs.scipy.org/doc/scipy/reference/_images/interpolate-3.png"

filename = url.split("/")[-1]
# print(filename)

datadir = os.path.join(os.path.dirname(__file__), "data", filename)

urllib.request.urlretrieve(url, datadir)

#########################################################################################
# urllib.request.urlretrieve("https://cofbkufrjwoilfgihmxyke.coursera-apps.org/view/week3/Car%20detection%20for%20Autonomous%20Driving/images/0012.jpg", "localname.jpg")
# https://cofbkufrjwoilfgihmxyke.coursera-apps.org/files/week3/Car%20detection%20for%20Autonomous%20Driving/images/0012.jpg


import requests

url1 = "https://cofbkufrjwoilfgihmxyke.coursera-apps.org/files/week3/Car%20detection%20for%20Autonomous%20Driving/images/0021.jpg"
r = requests.get(url, allow_redirects=False, stream=True)  # allow_redirects=True,
if r.status_code == 200:
    with open(
        os.path.join(os.path.dirname(__file__), "data", url1.split("/")[-1]), "wb"
    ) as f:
        f.write(r.content)
        # r.raw.decode_content = True
        # shutil.copyfileobj(r.raw, f)


url2 = "https://cofbkufrjwoilfgihmxyke.coursera-apps.org/files/week3/Car%20detection%20for%20Autonomous%20Driving/model_data/yolo.h5"
r = requests.get(url2, allow_redirects=False, stream=True)  # allow_redirects=True,
open(os.path.join(os.path.dirname(__file__), "data", url2.split("/")[-1]), "wb").write(
    r.content
)

# print(r.headers)
# print(r.headers.get('content-type'))
