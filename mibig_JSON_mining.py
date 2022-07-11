prod_activ = pd.read_csv('products_activity.csv', on_bad_lines='skip')
print(prod_activ)

BGCs_with_activity_list = []

for dirpath, dirnames, files in os.walk('../'):
    for file_name in files:
        if file_name.endswith('json'):
            with open(file_name,'r') as json_read:
                bgc_json = json.load(json_read)
                #print(bgc_json['cluster']['compounds'])
                comp = bgc_json['cluster']['compounds']
                for i in comp:
                    product = i['compound']
                    #print(product)
                    try:
                        const = prod_activ.loc[(prod_activ['cluster_name']).str.contains(str(product))]
                        if const.empty == False:
                            print(file_name)
                            print(const)
                            print('\n')
                            gbk_file_name = file_name.replace('.json','.gbk')
                            BGCs_with_activity_list.append(gbk_file_name)
                            print(gbk_file_name)
                            print('\n')
                        else:
                            print('empty_df')
                    except:
                        print("Error")
                #print('\n')
            #print(file_name)
print(BGCs_with_activity_list)
