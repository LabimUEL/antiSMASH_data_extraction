prod_activ = pd.read_csv('products_activity.csv', on_bad_lines='skip')
print(prod_activ)

for dirpath, dirnames, files in os.walk('./'):
    for file_name in files:
        if file_name.endswith('json'):
            with open(file_name,'r') as json_read:
                bgc_json = json.load(json_read)
                print(file_name)
                comp = bgc_json['cluster']['compounds']
                for i in comp:
                    product = i['compound']
                    print(product)
                print('\n')
