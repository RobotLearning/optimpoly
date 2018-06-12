import io
import json
import sys

file_json = sys.argv[1] #'../dmp.json'

with open(file_json,'r') as handle:
    parsed_json = json.load(handle)
    print json.dumps(parsed_json, indent = 4, sort_keys = True)

#parsed_json['tau'] = 1.0    

#with open(json_file, 'w') as json_file:
#    json.dump(json_decoded, json_file)

with io.open(file_json,'w',encoding='utf-8') as f:    
    f.write(json.dumps(parsed_json,indent=4,ensure_ascii=False,sort_keys=True))
    
