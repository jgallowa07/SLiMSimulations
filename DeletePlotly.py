import plotly
import plotly.plotly as py

import json
import requests
from requests.auth import HTTPBasicAuth

username = 'jgallowa' # Replace with YOUR USERNAME
api_key = '5migFKwWRcadU57uie1e' # Replace with YOUR API KEY

auth = HTTPBasicAuth(username, api_key)
headers = {'Plotly-Client-Platform': 'python'}

plotly.tools.set_credentials_file(username=username, api_key=api_key)
x = int(input("Which Graph?"))
numGraph = x

fid = username+':'+str(numGraph)

requests.post('https://api.plot.ly/v2/files/'+fid+'/trash', auth=auth, headers=headers)
