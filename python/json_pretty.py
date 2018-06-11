import io, json

file_json = '../dmp.json'

with open(file_json,'r') as handle:
    parsed_json = json.load(handle)
    print json.dumps(parsed_json, indent = 4, sort_keys = True)

with io.open(file_json,'w',encoding='utf-8') as f:
    f.write(json.dumps(parsed_json,indent=4,ensure_ascii=False,sort_keys=True))
    
